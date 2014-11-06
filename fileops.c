/*
 *  File operations routines
 */

#include "fileops.h"

/*****************************************************************************
* readims:  read image data from file with name "file"			     		*
*****************************************************************************/

int readims(char **filenames, maininfo* info, float *data){
    FILE *inputfile;
    size_t num_pixels = info->xsize * info->ysize;
    int i;

    for(i = 0; i < info->nrofframes; i++){
		inputfile = fopen(filenames[i], "r");

		if (inputfile != NULL) {
			/* Get offset out of the way in case it is present */
			fseek(inputfile, info->headeroffset, SEEK_SET);

			/* Read in the data in big chunks */
			fread(data + i*num_pixels, num_pixels * sizeof(float), 1, inputfile);

			fclose(inputfile);

			/* return success */
			return GO;
		} else {
			fprintf(stderr, "ERROR: file %s not found.\n", filenames[i]);
			return NOGO;
		}
    }

    return(GO);
}


/*****************************************************************************
* savefloat:  write float image data to file							     *
*****************************************************************************/

int savefloat(char filename[], long size, maininfo info, float *data)
{
    FILE *savefile = fopen(filename, "w");

    if (savefile != NULL) {
		/* Write data in one big chunks */
		fwrite(data, size * sizeof(float), 1, savefile);

		fclose(savefile);

		/* return success */
		return GO;
    } else {
    	fprintf(stderr, "ERROR: reconstruction could not be saved.\n");
		/* return NO success */
		return NOGO;
    }

    return(NOGO);
}
