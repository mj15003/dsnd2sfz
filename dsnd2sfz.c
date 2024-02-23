/**************************************************************************************/
/*    Utility for converting dsnd files of 2Box to WAV+sfz for use in Linuxsampler    */
/**************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>
#include <sndfile.h>
#include <errno.h>
#include <sys/stat.h>

#define MINUS_96DB 1.584893e-5f
#define HIT_LIMIT 255
#define ZONE_LIMIT 3

typedef enum {
	DSND_READ = 0,
	DSND_WAIT = 1
} DSND_WINDOW_STATE;

struct DSND_WINDOW {
	DSND_WINDOW_STATE state;
	sf_count_t start;
	sf_count_t len;
	float_t sum;
} win;

/* Struct for storing DSND zone info */

typedef struct {
	uint8_t	hit_n;				/* number of drum hits each will be separate WAV file*/
	sf_count_t start[HIT_LIMIT];		/* start index each hit */
	sf_count_t len[HIT_LIMIT];		/* length each hit in sndfile items */
	uint16_t num[HIT_LIMIT];		/* global hit identification number */
	uint8_t mvel[HIT_LIMIT];		/* calculated MIDI velocities */
	float_t peak_ampl[HIT_LIMIT];		/* maximum amplitude of each hit (dB)*/
	float_t max_peak_ampl;			/* maximum amplitude of zone as a whole (dB) */
	float_t min_peak_ampl;			/* minimum of peak amplitudes (dB)*/
	float_t dyn_range;			/* dynamic range of the zone in (dB) */
} DSND_ZONE_INFO;

/* Struct for storing DSND file info */

struct DSND_INFO {
	uint16_t hit_n;				/* number of drum hits each will be separate WAV file*/
	uint8_t zone_n;				/* number of separate zones (head, rim, sidestick) max. 3*/
	float_t peak_ampl;			/* overall maximum amplitude (dB)*/
	float_t min_peak_ampl;			/* minimum of peak amplitudes for a whole file (dB)*/
	float_t dyn_range;			/* overall dynamic range in dB between strongest and lightest hit amplitude */
	DSND_ZONE_INFO zones[ZONE_LIMIT];	/* array of zones in DSND file */
} dsnd_info;

SNDFILE *dsndf;			/* input DSND file */
float_t *samples;		/* input sample data */
SNDFILE *sf;			/* output WAV file(s) */

void signal_handler(int sig)
{
	fprintf(stderr, "\nsignal received, exiting ...\n");
	exit(EXIT_SUCCESS);
}

void show_usage()
{
	puts("\ndsnd2sfz -f dsnd-file-name [options]\n");
	puts("   -a dry run write no files");
	puts("   -b zone boundary volume difference (12 dB)");
	puts("   -d output directory");
	puts("   -h help");
	puts("   -n note (36)");
	puts("   -i in-threshold (60=-60 dB)");
	puts("   -o out-threshold (96=-96 dB)");
	puts("   -v number of velocity layers (5)");
	puts("   -w window length in samples (200)\n");
}

void cleanup()
{
	sf_close(dsndf);
	free(samples);
}

int main (int argc, char *argv[])
{
	/*******************/
	/* Initialization  */
	/*******************/

	char *dname;				/* Input DSND file name*/
	char *outdir;				/* Output directory */
	uint8_t note = 36;			/* MIDI note for SFZ file */
	uint8_t vel = 5;			/* number of requested velocity levels */
	float_t inthr = 1e-3f;			/* input threshold (-60 dB) */
	float_t outthr = MINUS_96DB;		/* output threshold (-96 dB) */
	float_t zbound = 12.0f;			/* zone boundary difference -> between peak amplitudes of the last hit in previous zone
						   and the first hit in the new zone */
	uint8_t dryrun = 0;			/* do analysis and write output, 1=analysis only */
	sf_count_t wlen = 200;			/* window length in frames (samples) */

	opterr = 0;
	int opt,err;
	char optok = 0;
	while ((opt = getopt(argc, argv, "ab:d:f:hn:i:o:v:w:")) != -1) {
		switch (opt) {
			case 'a':
				dryrun = 1;
				break;
			case 'b':
				zbound = strtof(optarg,NULL);
				break;
			case 'd':
				outdir = optarg;
				break;
			case 'f':
				dname = optarg;
				++optok;
				break;
			case 'h':
				show_usage();
				exit(EXIT_SUCCESS);
			case 'n':
				note = atoi(optarg);
				break;
			case 'i':
				inthr = powf(10.0f,-strtof(optarg, NULL)/20.0f);
				break;
			case 'o':
				outthr = powf(10.0f,-strtof(optarg, NULL)/20.0f);
				break;
			case 'v':
				vel = atoi(optarg);
				break;
			case 'w':
				wlen = atoi(optarg);
				break;
			default: /* '?' */
				show_usage();
				exit(EXIT_SUCCESS);
		}
	}

	if (argc < 3 || optok < 1) {
		show_usage();
		exit(EXIT_SUCCESS);
	}

	if (dname == NULL) {
		show_usage();
		exit(EXIT_FAILURE);
	}

	err = atexit(cleanup);
	if (err != 0)
	{
		fprintf(stderr, "cannot set exit function\n");
		exit(EXIT_FAILURE);
	}

	/************************************************************************/
	signal(SIGTERM, signal_handler);
	signal(SIGINT, signal_handler);
	/************************************************************************/
	/************************************************************************/
	/* Set parameters of sound files DSND,WAV                               */
	/************************************************************************/
	SF_INFO sf_info;	/* WAV file information */

	sf_count_t dsnd_items;	/* size of DSND sample data in items */
	size_t dsnd_bytes;	/* size of DSND sample data in bytes */

	sf_info.format = 0; /* other fields of SF_INFO will be filled in by library */

	if ((dsndf = sf_open(dname,SFM_READ,&sf_info)) == NULL) {
		char errstr[256];
		sf_error_str (0, errstr, sizeof (errstr) - 1);
		fprintf (stderr, "cannot open sndfile for input (%s)\n", errstr);
		exit (EXIT_FAILURE);
	} else {
		dsnd_bytes = sf_info.frames * sf_info.channels * sizeof(float_t);
		win.state = DSND_READ;
		win.sum = FLT_MIN;
		win.len = wlen * sf_info.channels;	/* length of the window is in items (frames * channels) */
		puts("DSND file parameters\n");
		printf("   name            : %s\n",dname);
		printf("   frames          : %u\n",sf_info.frames);
		printf("   samplerate      : %d\n",sf_info.samplerate);
		printf("   channels        : %d\n\n",sf_info.channels);
	}

	/* Allocate required memory and read DSND file*/
	samples = malloc(dsnd_bytes);
	dsnd_items = sf_read_float(dsndf, samples, sf_info.frames * sf_info.channels);
	printf("Read %u items from DSND file.\n\n",dsnd_items);

	/* Get overall DSND file statistics */

	dsnd_info.hit_n = 0;
	dsnd_info.zone_n = 1;
	dsnd_info.zones[0].peak_ampl[0] = -96.0f;

	/* Initialize structures */
	DSND_ZONE_INFO *zone = dsnd_info.zones;
	uint8_t hit = 0;

	/* Iterate through samples and gather stats and information */

	for (sf_count_t i=0;i<dsnd_items;i+=sf_info.channels) {
		/* we are evaluating channel 1 only */
		float_t s_abs_val = fabsf(samples[i]);
		float_t sdbv = 20 * log10f(s_abs_val);

		if (win.state == DSND_READ) {	/* reading sample data*/

			win.sum += sf_info.channels * s_abs_val / win.len; /* window length is in items we need frames */

			if (sdbv > zone->peak_ampl[hit]) zone->peak_ampl[hit] = sdbv;

			zone->len[hit] += sf_info.channels;	/* length is in items not frames*/

			if (zone->len[hit] >= win.len) {

				/* window contains 'len' samples and we can subtract as well*/
				if (zone->len[hit] == win.len) win.start = zone->start[hit];
				else win.start += sf_info.channels;

				win.sum -= sf_info.channels * fabsf(samples[win.start]) / win.len;

				if (win.sum < outthr) {
					/* End of the hit detected */
					win.state = DSND_WAIT;

					zone->num[hit] = ++dsnd_info.hit_n;
					zone->mvel[hit] = 64;

					if (hit > 0) {
						if (zone->peak_ampl[hit] > zone->peak_ampl[hit-1] + zbound) {
							/* Zone boundary detected */
							if (dsnd_info.zone_n < ZONE_LIMIT) {
								DSND_ZONE_INFO *nzone = zone+1;
								++dsnd_info.zone_n;
								/* Move last hit to the next zone*/
								nzone->start[0] = zone->start[hit];
								nzone->len[0] = zone->len[hit];
								nzone->num[0] = zone->num[hit];
								nzone->peak_ampl[0] = zone->peak_ampl[hit];
								nzone->mvel[0] = zone->mvel[hit];
								nzone->max_peak_ampl = zone->peak_ampl[hit];
								/* clear data for the hit in old zone*/
								zone->start[hit] = 0;
								zone->len[hit] = 0;
								zone->num[hit] = 0;
								zone->mvel[hit] = 64;
								zone->peak_ampl[hit] = -96.0f;

								zone = nzone;
								nzone->hit_n = 1;
								hit = 1;

								/* clear data for a new hit in the new zone*/
								zone->start[hit] = 0;
								zone->len[hit] = 0;
								zone->peak_ampl[hit] = -96.0f;
								zone->mvel[hit] = 64;

							} else {
								/* We already have ZONE_LIMIT zones -> no place to store data*/
								printf("Limit of %u zones has been reached.\n",ZONE_LIMIT);
								exit(EXIT_FAILURE);
							}
						} else {
						/* We are moving to the next hit in the same zone*/
							++hit;
							if (hit < HIT_LIMIT) {
								++zone->hit_n;
								zone->peak_ampl[hit] = -96.0f;
								zone->len[hit] = 0;
							} else {
								/* We already have HIT_LIMIT hits -> no place to store data*/
								printf("Limit of %u hits has been reached.\n",HIT_LIMIT);
								exit(EXIT_FAILURE);
							}
						}
					} else {
						/* We are moving to the next hit in the same zone*/
						++hit;
						++zone->hit_n;
						zone->peak_ampl[hit] = -96.0f;
						zone->len[hit] = 0;
					}
				}
			}
		} else {
			/* Start of the hit detected */
			if (s_abs_val > inthr) {
				win.state = DSND_READ;
				win.sum = sf_info.channels * s_abs_val / win.len;	/* window length is in items we need frames */
				zone->start[hit] = i;
			}
		}
	}

	/* Calculate dynamic range and MIDI velocities */
	dsnd_info.min_peak_ampl = 0.0f;
	dsnd_info.peak_ampl = -96.0f;		/* 16 Bit */

	for (uint8_t z=0; z < dsnd_info.zone_n; ++z) {
		zone = &dsnd_info.zones[z];
		zone->min_peak_ampl = 0.0f;
		zone->max_peak_ampl = -96.0f;

		for (uint8_t h=0; h < zone->hit_n; ++h) {
			if (zone->peak_ampl[h] < zone->min_peak_ampl) zone->min_peak_ampl = zone->peak_ampl[h];
			if (zone->peak_ampl[h] > zone->max_peak_ampl) zone->max_peak_ampl = zone->peak_ampl[h];
		}

		zone->dyn_range = zone->max_peak_ampl - zone->min_peak_ampl;

		for (uint8_t h=0; h < zone->hit_n; ++h) {
			zone->mvel[h] = (zone->peak_ampl[h] - zone->min_peak_ampl) * 127 / zone->dyn_range;
			if (zone->mvel[h] < 1) zone->mvel[h] = 1;
		}

		if (zone->max_peak_ampl > dsnd_info.peak_ampl) dsnd_info.peak_ampl = zone->max_peak_ampl;
		if (zone->min_peak_ampl < dsnd_info.min_peak_ampl) dsnd_info.min_peak_ampl = zone->min_peak_ampl;

	}

	dsnd_info.dyn_range = dsnd_info.peak_ampl - dsnd_info.min_peak_ampl;

	/* Print Zones and Hit table */

	puts("DSND file statistics\n");

	for (uint8_t z=0; z < dsnd_info.zone_n; ++z) {
		zone = &dsnd_info.zones[z];
		puts("+---+-----------+-----------+--------+----------+");
		puts("| N | Peak (dB) |   Start   | Length | MIDI vel |");
		puts("+---+-----------+-----------+--------+----------+");
		for (uint8_t h=0; h < zone->hit_n; ++h) {
			printf("|%3u|",zone->num[h]);
			printf("%11.2f|",zone->peak_ampl[h]);
			printf("%11u|",zone->start[h]/sf_info.channels);
			printf("%8u|",zone->len[h]/sf_info.channels);
			printf("%10u|\n",zone->mvel[h]);
		puts("+---+-----------+-----------+--------+----------+");
		}
	}
	printf("   Number of zones : %u\n",dsnd_info.zone_n);
	printf("   Number of hits  : %u\n\n",dsnd_info.hit_n);
	printf("   Peak amplitude  : %.2f dB\n",dsnd_info.peak_ampl);
	printf("   Dynamic range   : %.2f dB\n\n",dsnd_info.dyn_range);

	if (dryrun) exit(EXIT_SUCCESS);

	/************************************************************************/

	FILE *sfzfile;	/* Pointer to output SFZ file */

	char name_base[256];		/* name base for output file naming*/
	char wavname[256];		/* name for output WAV subor */
	char sfzname[256];		/* name for output SFZ file */

	if (chdir(outdir) < 0) {
		if (mkdir(outdir,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) < 0) {
			printf("Cannot make directory %s errno %d",outdir,errno);
		} else chdir(outdir);
	}

	/* Create the base name string from DSND file name */

	size_t dl = strlen(dname);
	uint8_t lslash = 255;
	uint8_t ldot = 0;
	for (uint8_t i = 0;i < dl;++i) {
		if (dname[i] == '/') lslash = i;
		else if (dname[i] == '.') ldot = i;
	}

	uint8_t n = 0;
	for (uint8_t i = lslash+1;i < ldot;++i) {
		name_base[n] = dname[i];
		++n;
	}
	name_base[n] = 0;
	/*****************************************************/

	sprintf(sfzname,"%s.sfz",name_base);
	printf("\nCreating SFZ file \'%s\'\n",sfzname);

	puts("\nWriting output WAV files.\n");

	/* Open SFZ file and write header */
	if ((sfzfile = fopen(sfzname,"w")) == NULL) {
		printf ("cannot open sfzfile for write (%s)\n", sfzname);
		exit (EXIT_FAILURE);
	}

	fputs("//***************************************************\n",sfzfile);
	fputs("//           generated by dsnd2sfz\n",sfzfile);
	fputs("//***************************************************\n",sfzfile);
	fputs("//***************************************************\n",sfzfile);
	fprintf(sfzfile,"//   Source file : %s\n",dname);
	fputs("//***************************************************\n",sfzfile);

	/***********************************************/
	/*  Generating SFZ file and exporting WAVs     */
	/***********************************************/
	fputs("<global> loop_mode=one_shot\n",sfzfile);

	for (uint8_t z=0; z < dsnd_info.zone_n; ++z) {
		zone = &dsnd_info.zones[z];

		uint8_t mvs = 128 / vel;	/* MIDI velocity step */
		int8_t hivel = 127;
		int8_t lovel = 128 - mvs;
		float_t gvmax = -96.0f;		/* group maximal volume */
		float_t gvmin = 0.0f;		/* group minimal volume */

		for (uint8_t v=0; v < vel; ++v) {

			uint8_t vlc = 0;	/* number of samples in selected velocity layer */

			/* iterating through hits and exporting respective WAVs */
			for (uint8_t h=0; h < zone->hit_n; ++h) {

				if (zone->mvel[h] <= hivel && zone->mvel[h] >= lovel) {

					++vlc;

					if (gvmax < zone->peak_ampl[h]) gvmax = zone->peak_ampl[h];
					if (gvmin > zone->peak_ampl[h]) gvmin = zone->peak_ampl[h];

					sprintf(wavname,"%s-%03u.wav",name_base,zone->num[h]);

					if ((sf = sf_open(wavname,SFM_WRITE,&sf_info)) == NULL) {
						char errstr[256];
						sf_error_str (0, errstr, sizeof (errstr) - 1);
						printf("cannot open sndfile for output (%s)\n", errstr);
						exit (EXIT_FAILURE);
					}

					/* Write the hit to output WAV file */
					sf_count_t owrt = sf_write_float(sf, &samples[zone->start[h]], zone->len[h]);
					sf_close(sf);
					printf("%s -> %u frames\n",wavname,owrt/sf_info.channels);
				}
			}

			if (vlc < 1) {
				lovel -= mvs;
				gvmax = -96.0f;
				gvmin = 0.0f;
				continue;
			}

			float_t mrs = 1.0f / vlc;	/* MIDI round-robin step */
			float_t lorand = 0.0f;
			float_t hirand = mrs;
			float_t vcr = powf(10.0f, (gvmin - gvmax)/20.0f);

			fprintf(sfzfile,"// Target volume = %.3f dB",gvmax);
			if (hivel == 127)
				fprintf(sfzfile,"\n<group> key=%u lovel=%u amp_velcurve_%u=%.3f\n",note,lovel,lovel,vcr);
			else if (lovel > 1)
				fprintf(sfzfile,"\n<group> key=%u lovel=%u hivel=%u amp_velcurve_%u=%.3f amp_velcurve_%u=1\n",note,lovel,hivel,lovel,vcr,hivel);
			else
				fprintf(sfzfile,"\n<group> key=%u hivel=%u amp_velcurve_1=%.3f amp_velcurve_%u=1\n",note,hivel,vcr,hivel);

			/* write regions in SFZ file for selected velocity layer group */
			uint8_t vlc2 = 0;	/* number of samples in selected veloctiy layer */

			for (uint8_t h=0; h < zone->hit_n; ++h) {

				if (zone->mvel[h] <= hivel && zone->mvel[h] >= lovel) {

					sprintf(wavname,"%s-%03u.wav",name_base,zone->num[h]);
					float_t gain = gvmax - zone->peak_ampl[h];

					if (vlc2 == 0) {
						fprintf(sfzfile,"<region> sample=%s volume=%.3f hirand=%.3f\n",wavname,gain,hirand);
					} else if (vlc2 < vlc-1) {
						fprintf(sfzfile,"<region> sample=%s volume=%.3f lorand=%.3f hirand=%.3f\n",wavname,gain,lorand,hirand);
					} else {
						fprintf(sfzfile,"<region> sample=%s volume=%.3f lorand=%.3f\n",wavname,gain,lorand);
					}

					++vlc2;
					lorand += mrs;
					hirand += mrs;
				}
			}

			hivel = lovel - 1;
			lovel = hivel - mvs;

			gvmax = -96.0f;
			gvmin = 0.0f;

			if (lovel < 1) lovel = 1;
		}
		++note;
	}

	fclose(sfzfile);
	printf("SFZ file \'%s\' has been created.\n",sfzname);

	/*******************************************/
	/*               Closing                   */
	/*******************************************/
	exit(EXIT_SUCCESS);
}
