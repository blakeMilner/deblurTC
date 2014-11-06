#ifndef OPTS_H
#define OPTS_H


/* default parameter values if no option specified */
#define DEFAULT_SAVEFILE		"TCaverage"
#define DEFAULT_XSIZE			0
#define DEFAULT_YSIZE			0
#define DEFAULT_HEADERSIZE		0
#define DEFAULT_SFS 			0.0
#define DEFAULT_BSAXISLENGTH 		20
#define DEFAULT_BSLENGTH 		20
#define DEFAULT_MAXRAD 			-1
#define DEFAULT_MAXITER 		10
#define DEFAULT_SNR 			1.0
#define DEFAULT_EXP 			0.2
#define DEFAULT_APODLIM			0.0

/* command line options for getopt() */
static struct option options[] = {
	 { "output",		required_argument,	NULL,	'o' },
	 { "xsize",		required_argument,	NULL,	'x' },
	 { "ysize",		required_argument,	NULL,	'y' },
	 { "headersize", 	required_argument,	NULL,	'h' },
	 { "subsize",		required_argument,	NULL,	's' },
	 { "bs1length",		required_argument,	NULL,	'v' },
	 { "bs2length",		required_argument,	NULL,	'w' },
	 { "maxrad",		required_argument,	NULL,	'r' },
	 { "maxiter",		required_argument,	NULL,	'p' },
	 { "snr",		required_argument,	NULL,	'n' },
	 { "weightexp"	,	required_argument,	NULL,	'e' },
	 { "apod",	 	required_argument,	NULL,	'a' }
};


#endif /* OPTS_H */
