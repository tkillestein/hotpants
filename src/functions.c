#include<stdio.h>
#include<string.h>
#include<math.h>
#include<stdlib.h>
#include<fitsio.h>
#include<ctype.h>

#include "defaults.h"
#include "globals.h"
#include "functions.h"

/*
  
Some of these subroutines appear originally in code created by Gary
Bernstein for psfmatch, but have been modified and/or rewritten for
the current software package.  In particular, the determination of
image noise statistics has been taken directly from the psfmatch
code.

08/20/01 acbecker@lucent.com

*/

/**
 * @brief Load an (X, Y) source position list from a whitespace-delimited text
 *        file into the global xcmp / ycmp arrays.
 *
 * @details Reads the first two columns of each non-alphabetic line as 1-indexed
 * (X, Y) coordinates and converts them to 0-indexed by subtracting 1.  If
 * @p cmpfileflag is non-zero the first line of the file is skipped (e.g. a
 * header line).  Reallocates xcmp and ycmp in chunks of 100 entries.
 *
 * @param filename     Path to the input text file.
 * @param cmpfileflag  If non-zero, skip the first line before parsing.
 */
/* armin */
void loadxyfile(char *filename, int cmpfileflag){
    FILE *xyfile;
    int Nalloc,c;
    char line[SCRLEN];
    xyfile   = fopen(filename, "r");
    fprintf(stderr, "WARNING : INPUT FORMAT HARDCODED : X=1, Y=2; 1-indexed coordinates\n");
    if (cmpfileflag) {
        for (;;){
            c=getc(xyfile);
            if (c== EOF) break;	 
            if (c=='\n') break;	  
        }
    }
    Nalloc=100;
    if ( !(xcmp = (float *)malloc(Nalloc*sizeof(float))) ||
         !(ycmp = (float *)malloc(Nalloc*sizeof(float))) )
        exit(1);
    Ncmp=0;
    while(fgets(line, SCRLEN, xyfile) != NULL) {
        /* reallocate if necessary */
        if (Ncmp==Nalloc){
            Nalloc+=100;
            xcmp  = (float *)realloc(xcmp, Nalloc *sizeof(float));
            ycmp  = (float *)realloc(ycmp, Nalloc *sizeof(float));
        }
        
        /* skip poorly formatted lines... */
        if (isalpha(line[0])) continue;
        if (sscanf(line, "%f %f", &xcmp[Ncmp], &ycmp[Ncmp]) == 2) {
            Ncmp++;
        }
    }
    fclose(xyfile);
    
    /* take care of 1-indexing, our arrays are 0-indexed */
    for (c = 0; c < Ncmp; c++) {
        xcmp[c] -= 1;
        ycmp[c] -= 1;
    }
}

/**
 * @brief Write used, all, and skipped substamp positions to DS9-compatible
 *        region files.
 *
 * @details Creates (or appends to, for regioncounter > 0) three region files:
 *   - xyfilename      : substamps selected for the kernel fit (green boxes)
 *   - xyfilename.all  : all candidate substamps (yellow boxes)
 *   - xyfilename.skipped : substamps that were tried but rejected (red boxes)
 * Coordinates are offset by (xmin, ymin) to convert from local region to global
 * image coordinates, and incremented by 1 to convert from 0-indexed to
 * 1-indexed FITS/DS9 convention.
 *
 * @param stamps         Array of stamps with substamp lists.
 * @param nStamps        Number of stamps.
 * @param xmin           X offset of the local region within the full image.
 * @param ymin           Y offset of the local region within the full image.
 * @param regioncounter  Pass 0 to create (overwrite) the files, >0 to append.
 */
/*armin*/
void savexy(stamp_struct *stamps, int nStamps, long xmin, long ymin, int regioncounter){
    int  sscnt,istamp;
    FILE *xyfileused,*xyfileall,*xyfileskipped;
    char xyfilenameall[1000];
    char xyfilenameskipped[1000];
    
    snprintf(xyfilenameall,     sizeof(xyfilenameall),     "%s.all",     xyfilename);
    snprintf(xyfilenameskipped, sizeof(xyfilenameskipped), "%s.skipped", xyfilename);
    
    if (regioncounter==0) {
        xyfileused = fopen(xyfilename, "w");
        xyfileall  = fopen(xyfilenameall, "w");
        xyfileskipped = fopen(xyfilenameskipped, "w");
    } else {
        xyfileused = fopen(xyfilename, "a");
        xyfileall  = fopen(xyfilenameall, "a");
        xyfileskipped = fopen(xyfilenameskipped, "a");
    }
    
    fprintf(xyfileall,    "#%4s %4s\n", "X", "Y");    
    fprintf(xyfileused,   "#%4s %4s\n", "X", "Y");    
    fprintf(xyfileskipped,"#%4s %4s\n", "X", "Y");
    
    /* note single pixel offset for ds9 display compared to 0-index array */
    for (istamp = 0; istamp < nStamps; istamp++) {
        for (sscnt = 0; sscnt < stamps[istamp].nss; sscnt++) {
            /*fprintf(xyfileall,        " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
            fprintf(xyfileall,        "image;box(%4ld,%4ld,%d,%d) # color=yellow\n",
                    stamps[istamp].xss[sscnt] + xmin + 1,
                    stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);    
            
            if (sscnt == stamps[istamp].sscnt)
                /*fprintf(xyfileused,    " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
                fprintf(xyfileused,    "image;box(%4ld,%4ld,%d,%d) # color=green\n",
                        stamps[istamp].xss[sscnt] + xmin + 1,
                        stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);
            
            if (sscnt < stamps[istamp].sscnt)
                /*fprintf(xyfileskipped, " %4ld %4ld\n", stamps[istamp].xss[sscnt] + xmin, stamps[istamp].yss[sscnt] + ymin);*/
                fprintf(xyfileskipped, "image;box(%4ld,%4ld,%d,%d) # color=red\n",
                        stamps[istamp].xss[sscnt] + xmin + 1,
                        stamps[istamp].yss[sscnt] + ymin + 1, fwKSStamp, fwKSStamp);    
        }
    }
    
    fclose(xyfileskipped);
    fclose(xyfileall);
    fclose(xyfileused);
}

/**
 * @brief Allocate all dynamic sub-arrays within an array of stamp_struct
 *        objects.
 *
 * @details For each stamp, allocates:
 *   - stamps[i].vectors : (nCompKer + nBGVectors) × (fwKSStamp²) double arrays
 *     for the convolved kernel basis images and background basis images.
 *   - stamps[i].krefArea: fwKSStamp² double array for the reference-image patch.
 *   - stamps[i].mat     : nC × nC double array for the per-stamp normal-equation
 *     matrix.
 *   - stamps[i].xss, stamps[i].yss: integer arrays of length nKSStamps for
 *     substamp coordinates.
 *   - stamps[i].scprod  : nC double array for the right-hand-side vector.
 * All scalar fields are zeroed.
 *
 * @param stamps   Pre-allocated array of nStamps stamp_struct objects (the
 *                 outer array must already exist).
 * @param nStamps  Number of stamps to initialise.
 * @return 0 on success, 1 if any allocation fails.
 */
int allocateStamps(stamp_struct *stamps, int nStamps) {
    int stampIdx;
    int nbgVectors;

    nbgVectors = ((bgOrder + 1) * (bgOrder + 2)) / 2;

    if (stamps) {
        for (stampIdx = 0; stampIdx < nStamps; stampIdx++) {
            /* Allocate contiguous vector array: much more efficient than nested malloc */
            stamps[stampIdx].vectors = alloc_vector_array(nCompKer + nbgVectors, fwKSStamp * fwKSStamp);

            /* Reference area for kernel fitting */
            stamps[stampIdx].krefArea = (double *)xcalloc((size_t)fwKSStamp * fwKSStamp, sizeof(double));

            /* Fitting matrix: contiguous 2D array for better cache locality */
            stamps[stampIdx].mat = alloc_matrix_contiguous(nC, nC);

            /* Substamp center coordinates */
            stamps[stampIdx].xss = (int *)xcalloc(nKSStamps, sizeof(int));
            stamps[stampIdx].yss = (int *)xcalloc(nKSStamps, sizeof(int));

            /* Scalar products (right-hand side of normal equations) */
            stamps[stampIdx].scprod = (double *)xcalloc(nC, sizeof(double));

            /* Initialize stamp metadata */
            stamps[stampIdx].x0 = stamps[stampIdx].y0 = stamps[stampIdx].x = stamps[stampIdx].y = 0;
            stamps[stampIdx].nss = stamps[stampIdx].sscnt = 0;
            stamps[stampIdx].nx = stamps[stampIdx].ny = 0;
            stamps[stampIdx].sum = stamps[stampIdx].mean = stamps[stampIdx].median = 0;
            stamps[stampIdx].mode = stamps[stampIdx].sd = stamps[stampIdx].fwhm = 0;
            stamps[stampIdx].lfwhm = stamps[stampIdx].chi2 = 0;
            stamps[stampIdx].norm = stamps[stampIdx].diff = 0;
        }
    }
    return 0;
}


/**
 * @brief Initialise stamp statistics and discover or assign substamp centres
 *        within a stamp region.
 *
 * @details For each of the template and image stamp arrays (ctStamps, ciStamps),
 * if no substamps have been found yet, calls cutStamp() to extract the stamp
 * pixel data and getStampStats3() to compute sky background statistics (mode,
 * FWHM used as sky sigma).
 *
 * Substamp centres are then assigned in one of two modes:
 *  - getCenters == 1: calls getPsfCenters() to automatically locate bright,
 *    isolated PSF-like stars within the stamp, up to nKSStamps substamps.
 *  - getCenters == 0: assigns a single substamp at the stamp centre (or at a
 *    hard-coded position if hardX/hardY are non-zero) after validating it with
 *    checkPsfCenter().  The selected region is masked in mRData with
 *    FLAG_T_SKIP / FLAG_I_SKIP to prevent overlap with other stamps.
 *
 * @param sXMin      Left pixel bound of the stamp region in global image coords.
 * @param sXMax      Right pixel bound.
 * @param sYMin      Bottom pixel bound.
 * @param sYMax      Top pixel bound.
 * @param niS        In/out: index of the current image stamp in ciStamps.
 * @param ntS        In/out: index of the current template stamp in ctStamps.
 * @param getCenters 1 to run automatic PSF-centre search; 0 for single-centre
 *                   assignment.
 * @param rXBMin     X origin of the region buffer (offset from full image).
 * @param rYBMin     Y origin of the region buffer.
 * @param ciStamps   Array of image stamps being built.
 * @param ctStamps   Array of template stamps being built.
 * @param iRData     Full-frame image pixel data.
 * @param tRData     Full-frame template pixel data.
 * @param hardX      If non-zero, force the substamp X centre to this value.
 * @param hardY      If non-zero, force the substamp Y centre to this value.
 */
void buildStamps(int sXMin, int sXMax, int sYMin, int sYMax, int *niS, int *ntS,
                 int getCenters, int rXBMin, int rYBMin, stamp_struct *ciStamps, stamp_struct *ctStamps,
                 float *iRData, float *tRData, float hardX, float hardY) {
    
    int sPixX, sPixY;
    int xmax, ymax, k, l, nr2;
    int nss;
    /*int bbitt1=0x100, bbitt2=0x200, bbiti1=0x400, bbiti2=0x800;*/
    int bbitt1=FLAG_T_BAD, bbitt2=FLAG_T_SKIP, bbiti1=FLAG_I_BAD, bbiti2=FLAG_I_SKIP;
    float *refArea=NULL;
    double check;
    
    if (verbose >= 1)
        fprintf(stderr, "    Stamp in region : %d:%d,%d:%d\n",
                sXMin, sXMax, sYMin, sYMax);
    
    /* global vars */
    sPixX  = sXMax - sXMin + 1;
    sPixY  = sYMax - sYMin + 1;
    
    if (!(strncmp(forceConvolve, "i", 1)==0)) {
        
        if (ctStamps[*ntS].nss == 0) {
            refArea = (float *)calloc(sPixX*sPixY, sizeof(float));
            
            /* temp : store the whole stamp in refArea */
            cutStamp(tRData, refArea, rPixX, sXMin-rXBMin, sYMin-rYBMin,
                     sXMax-rXBMin, sYMax-rYBMin, &ctStamps[*ntS]);
            
            if ( !( getStampStats3(refArea, ctStamps[*ntS].x0, ctStamps[*ntS].y0, sPixX, sPixY,
                                   &ctStamps[*ntS].sum, &ctStamps[*ntS].mean, &ctStamps[*ntS].median,
                                   &ctStamps[*ntS].mode, &ctStamps[*ntS].sd, &ctStamps[*ntS].fwhm,
                                   &ctStamps[*ntS].lfwhm, 0x0, 0xffff, 3))) {
                /* pointless */
                if (verbose >= 1)
                    fprintf(stderr, "    Tmpl  xs : %4i ys : %4i  (sky,dsky = %.1f,%.1f)\n",
                            ctStamps[*ntS].x, ctStamps[*ntS].y, ctStamps[*ntS].mode, ctStamps[*ntS].fwhm);
            }
            free(refArea);
        }
    }
    
    if (!(strncmp(forceConvolve, "t", 1)==0)) {
        
        if (ciStamps[*niS].nss == 0) {
            refArea = (float *)calloc(sPixX*sPixY, sizeof(float));
            
            /* store the whole of the stamp in .refArea */
            cutStamp(iRData, refArea, rPixX, sXMin-rXBMin, sYMin-rYBMin,
                     sXMax-rXBMin, sYMax-rYBMin, &ciStamps[*niS]);
            
            if ( !( getStampStats3(refArea, ciStamps[*niS].x0, ciStamps[*niS].y0, sPixX, sPixY,
                                   &ciStamps[*niS].sum, &ciStamps[*niS].mean, &ciStamps[*niS].median,
                                   &ciStamps[*niS].mode, &ciStamps[*niS].sd, &ciStamps[*niS].fwhm,
                                   &ciStamps[*niS].lfwhm, 0x0, 0xffff, 3))) {
                
                /* pointless */
                /* buildSigMask(&ciStamps[*niS], sPixX, sPixY, misRData); */
                if (verbose >= 1)
                    fprintf(stderr, "    Image xs : %4i ys : %4i  (sky,dsky = %.1f,%.1f)\n",
                            ciStamps[*niS].x, ciStamps[*niS].y, ciStamps[*niS].mode, ciStamps[*niS].fwhm);
            }
            free(refArea);
        }
    }
    
    if (!(strncmp(forceConvolve, "i", 1)==0)) {
        nss = ctStamps[*ntS].nss;
        if (getCenters) {
            /* get potential centers for the kernel fit */
            getPsfCenters(&ctStamps[*ntS], tRData, sPixX, sPixY, tUKThresh, bbitt1, bbitt2);
            if (verbose >= 1)
                fprintf(stderr, "    Tmpl     : scnt = %2i nss = %2i\n",
                        ctStamps[*ntS].sscnt, ctStamps[*ntS].nss);
            
        } else {
            /* don't increment ntS inside subroutine, BUT MAKE SURE YOU DO IT OUTSIDE! */
            if (nss < nKSStamps) {
                
                if (hardX)
                    xmax = (int)(hardX);
                else
                    xmax = sXMin + (int)(fwStamp/2);
                
                if (hardY)
                    ymax = (int)(hardY);
                else
                    ymax = sYMin + (int)(fwStamp/2);
                
                check = checkPsfCenter(tRData, xmax-ctStamps[*ntS].x0, ymax-ctStamps[*ntS].y0, sPixX, sPixY,
                                       ctStamps[*ntS].x0, ctStamps[*ntS].y0, tUKThresh,
                                       ctStamps[*ntS].mode, 1. / ctStamps[*ntS].fwhm,
                                       0, 0, bbitt1 | bbitt2 | 0xbf, bbitt1);
                
                /* its ok? */
                if (check != 0.) {
                    /* globally mask out the region around this guy */
                    for (l = ymax-hwKSStamp; l <= ymax+hwKSStamp; l++) {
                        /*yr2 = l + ctStamps[*ntS].y0;*/
                        
                        for (k = xmax-hwKSStamp; k <= xmax+hwKSStamp; k++) {
                            /*xr2 = k + ctStamps[*ntS].x0;*/
                            nr2 = l+rPixX*k;
                            
                            /*if ((k > 0) && (k < sPixX) && (l > 0) && (l < sPixY))*/
                            if (nr2 >= 0 && nr2 < rPixX*rPixY)
                                mRData[nr2] |= bbitt2;
                        }
                    }
                    
                    ctStamps[*ntS].xss[nss] = xmax;
                    ctStamps[*ntS].yss[nss] = ymax;	    
                    ctStamps[*ntS].nss += 1;
                    if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d\n", nss, xmax, ymax);
                }
            }
        }
    }
    
    if (!(strncmp(forceConvolve, "t", 1)==0)) {
        nss = ciStamps[*niS].nss;
        
        if (getCenters) {
            /* get potential centers for the kernel fit */
            getPsfCenters(&ciStamps[*niS], iRData, sPixX, sPixY, iUKThresh, bbiti1, bbiti2);
            if (verbose >= 1)
                fprintf(stderr, "    Image    : scnt = %2i nss = %2i\n",
                        ciStamps[*niS].sscnt, ciStamps[*niS].nss);
        } else {
            /* don't increment niS inside subroutine, BUT MAKE SURE YOU DO IT OUTSIDE! */
            if (nss < nKSStamps) {
                if (hardX)
                    xmax = (int)(hardX);
                else
                    xmax = sXMin + (int)(fwStamp/2);
                
                if (hardY)
                    ymax = (int)(hardY);
                else
                    ymax = sYMin + (int)(fwStamp/2);
                
                check = checkPsfCenter(iRData, xmax-ciStamps[*niS].x0, ymax-ciStamps[*niS].y0, sPixX, sPixY,
                                       ciStamps[*niS].x0, ciStamps[*niS].y0, iUKThresh,
                                       ciStamps[*niS].mode, 1. / ciStamps[*niS].fwhm,
                                       0, 0, bbiti1 | bbiti2 | 0xbf, bbiti1);
                
                /* its ok? */
                if (check != 0.) {
                    /* globally mask out the region around this guy */
                    for (l = ymax-hwKSStamp; l <= ymax+hwKSStamp; l++) {
                        /*yr2 = l + ciStamps[*niS].y0;*/
                        
                        for (k = xmax-hwKSStamp; k <= xmax+hwKSStamp; k++) {
                            /*xr2 = k + ciStamps[*niS].x0;*/
                            /*nr2 = xr2+rPixX*yr2;*/
                            nr2 = l+rPixX*k;
                            
                            /*if ((k > 0) && (k < sPixX) && (l > 0) && (l < sPixY))*/
                            if (nr2 >= 0 && nr2 < rPixX*rPixY)
                                mRData[nr2] |= bbiti2;
                        }
                    }
                    
                    ciStamps[*niS].xss[nss] = xmax;
                    ciStamps[*niS].yss[nss] = ymax;	    
                    ciStamps[*niS].nss += 1;
                    if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d\n", nss, xmax, ymax);
                    
                }
            }
        }
    }
    
    return;
}

/**
 * @brief Copy a rectangular sub-region of a full-frame image into a contiguous
 *        buffer and record the stamp's position metadata.
 *
 * @details Copies pixels data[xMin..xMax][yMin..yMax] (inclusive) into the
 * pre-allocated array @p refArea, reindexed to [0..sxLen-1][0..syLen-1].  Also
 * sets stamp->x0, stamp->y0 to the region origin and stamp->x, stamp->y to the
 * region centre.
 *
 * @param data    Full-frame input image.
 * @param refArea Output buffer of size (xMax-xMin+1)*(yMax-yMin+1), pre-allocated
 *                by the caller.
 * @param dxLen   Row stride of @p data (i.e. full-image width, rPixX).
 * @param xMin    Left column index (inclusive) in @p data.
 * @param yMin    Bottom row index (inclusive) in @p data.
 * @param xMax    Right column index (inclusive) in @p data.
 * @param yMax    Top row index (inclusive) in @p data.
 * @param stamp   Stamp whose origin and centre fields are set on return.
 */
void cutStamp(float *data, float *refArea, int dxLen, int xMin, int yMin,
              int xMax, int yMax, stamp_struct *stamp) {
    
    int srcCol, srcRow;
    int x, y;
    int sxLen;

    sxLen = xMax - xMin + 1;
    /* NOTE, use srcRow <= yMax here, not srcRow < yMax, include yMax point */
    for (srcRow = yMin; srcRow <= yMax; srcRow++) {
        y = srcRow - yMin;

        for (srcCol = xMin; srcCol <= xMax; srcCol++) {
            x = srcCol - xMin;

            refArea[x+y*sxLen] = data[srcCol+srcRow*dxLen];
            /*fprintf(stderr, "%d %d %d %d %f\n", x, y, srcCol, srcRow, data[srcCol+srcRow*dxLen]); */
        }
    }
    stamp->x0 = xMin;
    stamp->y0 = yMin;
    stamp->x  = xMin + (xMax - xMin) / 2;
    stamp->y  = yMin + (yMax - yMin) / 2;
    
    return;
}

/**
 * @brief Extract the reference-image pixel patch for the current substamp into
 *        stamp->krefArea.
 *
 * @details The substamp centre (stamp->xss[sscnt], stamp->yss[sscnt]) is given
 * in global image coordinates.  This function converts to local stamp
 * coordinates using stamp->x0, stamp->y0, then copies the
 * fwKSStamp×fwKSStamp pixel region from iData into stamp->krefArea, which is
 * used later as the fitting target in build_scprod0() and getStampSig().
 * Masked (FLAG_INPUT_ISBAD) pixels contribute 0 to stamp->sum.
 * Returns 1 if all substamps are exhausted (sscnt >= nss), causing the stamp
 * to be rejected.
 *
 * @param stamp  Stamp being processed; stamp->sscnt selects the active substamp.
 * @param iData  Full-frame reference image in global pixel coordinates.
 * @return 0 on success, 1 if the stamp is out of valid substamps.
 */
int cutSStamp(stamp_struct *stamp, float *iData) {

    int substampCol, substampRow, k, nss, sscnt;
    int x, y, dy, xStamp, yStamp;
    float dpt;
    double sum = 0.;
    
    /*stamp->krefArea = (double *)realloc(stamp->krefArea, fwKSStamp*fwKSStamp*sizeof(double));*/
    dfset(stamp->krefArea, fillVal, fwKSStamp, fwKSStamp);
    
    nss    = stamp->nss;
    sscnt  = stamp->sscnt;
    xStamp = stamp->xss[sscnt] - stamp->x0;
    yStamp = stamp->yss[sscnt] - stamp->y0;
    
    /* have gone through all the good substamps, reject this stamp */
    if (sscnt >= nss) {
        if (verbose >= 2)
            fprintf(stderr, "    xs : %4i ys : %4i sig: %6.3f sscnt: %4i nss: %4i ******** REJECT stamp (out of substamps)\n",
                    stamp->x, stamp->y, stamp->chi2, sscnt, nss);
        else
            if (verbose >= 1)
                fprintf(stderr, "        Reject stamp\n");
        return 1;
    }
    /*
      fprintf(stderr, "    xs : %4i ys : %4i substamp (sscnt=%d, nss=%d) xss: %4i yss: %4i\n",
      stamp->x, stamp->y, sscnt, nss, stamp->xss[sscnt], stamp->yss[sscnt]);
    */
    if (verbose >= 1)
        fprintf(stderr, "    xss : %4i yss : %4i\n", stamp->xss[sscnt], stamp->yss[sscnt]);
    
    for (substampRow = yStamp - hwKSStamp; substampRow <= yStamp + hwKSStamp; substampRow++) {
        y  = substampRow - (yStamp - hwKSStamp);
        dy = substampRow + stamp->y0;

        for (substampCol = xStamp - hwKSStamp; substampCol <= xStamp + hwKSStamp; substampCol++) {
            x = substampCol - (xStamp - hwKSStamp);

            k   = substampCol+stamp->x0 + rPixX*dy;
            dpt = iData[k];

            stamp->krefArea[x+y*fwKSStamp] = dpt;
            sum += (mRData[k] & FLAG_INPUT_ISBAD) ? 0 : fabs(dpt);
        }
    }
    
    stamp->sum = sum;
    return 0;
}

/**
 * @brief Validate a candidate substamp centre and return a ranking score.
 *
 * @details Checks the fwKSStamp×fwKSStamp box centred on (imax, jmax) in the
 * stamp coordinate system for any disqualifying conditions:
 *  - Any pixel masked with @p bbit (overlap, saturation, etc.) → score = 0.
 *  - Any pixel exceeding @p hiThresh (saturation) → flagged with @p bbit1,
 *    score = 0.
 * If the box passes, the score is the sum of pixel values exceeding
 * kerFitThresh standard deviations above the sky background, providing a
 * brightness-based ranking so that brighter, more isolated PSF candidates are
 * preferred.
 *
 * @param iData    Full-frame (or stamp-region) image.
 * @param imax     Candidate centre x in stamp coordinates.
 * @param jmax     Candidate centre y in stamp coordinates.
 * @param xLen     Width of the stamp region.
 * @param yLen     Height of the stamp region.
 * @param sx0      X offset of the stamp region within the full image.
 * @param sy0      Y offset of the stamp region within the full image.
 * @param hiThresh Pixel value above which a star is considered saturated.
 * @param sky      Sky background level (mode of stamp pixel distribution).
 * @param invdsky  1 / sky_sigma (reciprocal of the sky noise, i.e. 1/FWHM).
 * @param xbuffer  Minimum x margin from stamp edge (typically 0).
 * @param ybuffer  Minimum y margin from stamp edge (typically 0).
 * @param bbit     Mask bits that disqualify a candidate (any overlap/bad flag).
 * @param bbit1    Mask bit to set on saturated pixels.
 * @return Brightness score (sum of above-threshold pixel values) if valid, 0 if
 *         the candidate is disqualified.
 */
double checkPsfCenter(float *iData, int imax, int jmax, int xLen, int yLen,
                      int sx0, int sy0,
                      double hiThresh, float sky, float invdsky,
                      int xbuffer, int ybuffer, int bbit, int bbit1) {

    /* note the imax and jmax are in the stamp coordinate system */

    int brk, checkRowIdx, checkColIdx;
    int xr2, yr2, nr2;
    double dmax2, dpt2;
    
    /* (5). after centroiding, re-check for fill/sat values and overlap (a little inefficient) */
    /*      zero tolerance for bad pixels! */
    brk   = 0;
    
    /* since we have zero tolerance for bad pixels, the sum
       of all the pixels in this box can be compared to the
       sum of all the pixels in the other boxes.  rank on
       this!  well, sum of high sigma pixels anyways...*/
    dmax2 = 0.;

    for (checkRowIdx = jmax-hwKSStamp; checkRowIdx <= jmax+hwKSStamp; checkRowIdx++) {
        if ((checkRowIdx < ybuffer) || (checkRowIdx >= yLen-ybuffer))
            continue; /* continue checkRowIdx loop */

        yr2 = checkRowIdx + sy0;

        for (checkColIdx = imax-hwKSStamp; checkColIdx <= imax+hwKSStamp; checkColIdx++) {
            if ((checkColIdx < xbuffer) || (checkColIdx >= xLen-xbuffer))
                continue; /* continue checkColIdx loop */

            xr2 = checkColIdx + sx0;
            nr2 = xr2+rPixX*yr2;

            if (mRData[nr2] & bbit) {
                brk = 1;
                dmax2 = 0;
                break; /* exit checkColIdx loop */
            }

            dpt2   = iData[nr2];

            if (dpt2 >= hiThresh) {
                mRData[nr2] |= bbit1;
                brk = 1;
                dmax2 = 0;
                break; /* exit checkColIdx loop */
            }

            if (( (dpt2 - sky) * invdsky) > kerFitThresh)
                dmax2 += dpt2;
        }
        if (brk == 1)
            break; /* exit checkRowIdx loop */
    }
    
    return dmax2;
}

/**
 * @brief Automatically locate up to nKSStamps bright, isolated PSF-like stars
 *        within a stamp region for use as kernel-fitting substamps.
 *
 * @details Iteratively lowers a detection threshold from hiThresh × dfrac down
 * to a floor of sky + kerFitThresh·sigma in steps of dfrac -= 0.2.  At each
 * threshold level, scans all unmasked pixels above the current threshold,
 * re-centres each candidate on the brightest pixel within an hwKSStamp box,
 * then calls checkPsfCenter() to validate the fwKSStamp region (zero tolerance
 * for bad pixels).  Valid candidates are scored by their above-threshold flux
 * sum and recorded in xloc/yloc/peaks arrays.  Accepted regions are masked with
 * bbit2 (FLAG_T_SKIP or FLAG_I_SKIP) to prevent overlap.
 *
 * After scanning, candidates are sorted by brightness using quick_sort() and
 * the top nKSStamps are stored in stamp->xss[], stamp->yss[].
 *
 * @param stamp     Stamp for which substamp centres are to be found; stamp->nss
 *                  is incremented for each accepted substamp.
 * @param iData     Full-frame (or region) image.
 * @param xLen      Width of the stamp region in pixels.
 * @param yLen      Height of the stamp region in pixels.
 * @param hiThresh  Saturation / upper threshold for valid pixels.
 * @param bbit1     Mask bit for saturated pixels (e.g. FLAG_T_BAD).
 * @param bbit2     Mask bit for accepted-substamp exclusion zones
 *                  (e.g. FLAG_T_SKIP).
 * @return 0 on success, 1 if no valid substamp could be found.
 */
int getPsfCenters(stamp_struct *stamp, float *iData, int xLen, int yLen, double hiThresh, int bbit1, int bbit2) {

    int stampRegionCol, stampRegionRow, centroidCol, centroidRow, nr, nr2, imax, jmax, xr, yr, xr2, yr2, sy0, sx0, xbuffer, ybuffer, i, j;
    double dmax, dmax2, dpt, dpt2, loPsf, floor;
    int *xloc, *yloc, *qs, pcnt, brk, bbit, fcnt;
    double *peaks;
    double sky, invdsky;
    
    float dfrac = 0.9;
    
    if (stamp->nss >= nKSStamps) {
        fprintf(stderr,"    no need for automatic substamp search...\n");
        return(0);
    }
    
    /* 0xff, but 0x40 is OK = 0xbf */
    bbit    = bbit1 | bbit2 | 0xbf;
    sky     = stamp->mode;
    invdsky = 1. / stamp->fwhm;
    
    /* the highest peak in the image, period, in case all else fails */
    /* default to center of stamp */
    sx0  = stamp->x0;
    sy0  = stamp->y0;
    
    /*
      we can go with no buffer here, but, we need to allocate extra size in the
      stamp->refArea to allow this to happen.  too much room.  removed
      refArea from stamp structure, just use actual array.
    */
    xbuffer = ybuffer = 0;
    
    /* as low as we will go */
    /* was a 2 here before version 4.1.6 */
    floor = sky + kerFitThresh * stamp->fwhm;
    
    qs    = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    xloc  = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    yloc  = (int *)calloc(xLen*yLen/(hwKSStamp), sizeof(int));
    peaks = (double *)calloc(xLen*yLen/(hwKSStamp), sizeof(double));
    
    brk  = 0;
    pcnt = 0;
    fcnt = 2 * nKSStamps; /* maximum number of stamps to find */
    
    while (pcnt < fcnt) {
        
        /* NOTE : we keep the previous iteration's matches if they were good */
        
        loPsf = sky + (hiThresh - sky) * dfrac;
        loPsf = (loPsf > floor) ? loPsf : floor;   /* last ditch attempt to get some usable pixels */
        
        /* (0). ignore near the edges */
        for (stampRegionRow = ybuffer; stampRegionRow < yLen-ybuffer; stampRegionRow++) {
            yr  = stampRegionRow + sy0;

            for (stampRegionCol = xbuffer; stampRegionCol < xLen-xbuffer; stampRegionCol++) {
                xr = stampRegionCol + sx0;
                nr = xr+rPixX*yr;
                
                /* (1). pixel already included in another stamp, or exceeds hiThresh */
                if (mRData[nr] & bbit)
                    continue;
                
                dpt = iData[nr];
                
                /* (2a). skip masked out value */
                if (dpt >= hiThresh) {
                    mRData[nr] |= bbit1;
                    continue;
                }
                
                /* (2b). don't want low sigma point here, sometimes thats all that passes first test */
                if (( (dpt - sky) * invdsky) < kerFitThresh)
                    continue;
                
                /* finally, a good candidate! */
                if (dpt > loPsf) {

                    dmax  = dpt;
                    imax  = stampRegionCol;
                    jmax  = stampRegionRow;
                    
                    /* center this candidate on the peak flux */
                    for (centroidRow = stampRegionRow-hwKSStamp; centroidRow <= stampRegionRow+hwKSStamp; centroidRow++) {
                        yr2 = centroidRow + sy0;

                        if ((centroidRow < ybuffer) || (centroidRow >= yLen-ybuffer))
                            continue; /* continue centroidRow loop */

                        for (centroidCol = stampRegionCol-hwKSStamp; centroidCol <= stampRegionCol+hwKSStamp; centroidCol++) {
                            xr2 = centroidCol + sx0;
                            nr2 = xr2+rPixX*yr2;

                            if ((centroidCol < xbuffer) || (centroidCol >= xLen-xbuffer))
                                continue; /* continue centroidCol loop */
                            
                            /* (3). no fill/sat values within checked region, and not overlapping other stamp */
                            if (mRData[nr2] & bbit)
                                continue; /* continue k loop */
                            
                            dpt2 = iData[nr2];
                            
                            /* (4a). no fill/sat values within checked region, and not overlapping other stamp */
                            if (dpt2 >= hiThresh) {
                                mRData[nr2] |= bbit1;
                                continue; /* continue k loop */
                            }
                            
                            /* (4b). not as strong a problem, just don't want low sigma peak */
                            if (( (dpt2 - sky) * invdsky) < kerFitThresh)
                                continue;
                            
                            /* record position and amp of good point centroid in stampRegionCol,stampRegionRow,dmax */
                            if (dpt2 > dmax) {
                                dmax = dpt2;
                                imax = centroidCol;
                                jmax = centroidRow;
                            }
                        }
                    }
                    
                    /* (5). after centroiding, re-check for fill/sat values and overlap (a little inefficient) */
                    /*      zero tolerance for bad pixels! */
                    
                    /* since we have zero tolerance for bad pixels, the sum
                       of all the pixels in this box can be compared to the
                       sum of all the pixels in the other boxes.  rank on
                       this!  well, sum of high sigma pixels anyways...*/
                    
                    dmax2 = checkPsfCenter(iData, imax, jmax, xLen, yLen, sx0, sy0,
                                           hiThresh, sky, invdsky, xbuffer, ybuffer, bbit, bbit1);
                    
                    if (dmax2 == 0.)
                        continue; /* continue stampRegionCol loop */
                    
                    /* made it! - a valid peak */
                    xloc[pcnt]    = imax;
                    yloc[pcnt]    = jmax;
                    peaks[pcnt++] = dmax2;
                    
                    /* globally mask out the region around this guy */
                    for (centroidRow = jmax-hwKSStamp; centroidRow <= jmax+hwKSStamp; centroidRow++) {
                        yr2 = centroidRow + sy0;

                        for (centroidCol = imax-hwKSStamp; centroidCol <= imax+hwKSStamp; centroidCol++) {
                            xr2 = centroidCol + sx0;
                            nr2 = xr2+rPixX*yr2;

                            if ((centroidCol > 0) && (centroidCol < xLen) && (centroidRow > 0) && (centroidRow < yLen))
                                mRData[nr2] |= bbit2;
                        }
                    }
                    
                    /* found enough stamps, get out! */
                    if (pcnt >= fcnt)
                        brk = 2;
                }
                if (brk == 2)
                    break; /* exit x loop */
            }
            if (brk == 2)
                break; /* exit y loop */
        }
        
        if (loPsf == floor)
            break;        /* have hit the floor, get out  */
        dfrac -= 0.2;
        
    }
    
    if (pcnt+stamp->nss < nKSStamps) {
        if (verbose >= 2) fprintf(stderr, "    ...only found %d good substamps by autosearch\n", pcnt);
    }
    else {
        if (verbose >= 2) fprintf(stderr, "    ...found %d good substamps by autosearch\n", pcnt);
    }
    
    /* coords are in region's pixels */
    
    /* no good full stamps */
    if (pcnt == 0) {
        /* changed in v4.1.6, don't accept center pixel, could be bad, worse than not having one! */
        if (verbose >= 2) fprintf(stderr, "    NO good pixels, skipping...\n");
        
        free(qs);
        free(xloc);
        free(yloc);
        free(peaks);
        return 1;
    }
    else {
        quick_sort (peaks, qs, pcnt);
        if (verbose >= 2) fprintf(stderr, "    Adding %d substamps found by autosearch\n", imin(pcnt, nKSStamps-stamp->nss));
        for (i = stamp->nss, j = 0; j < pcnt && i < nKSStamps; i++,j++) {
            stamp->xss[i] = xloc[qs[pcnt-j-1]] + sx0;
            stamp->yss[i] = yloc[qs[pcnt-j-1]] + sy0;
            if (verbose >= 2) fprintf(stderr,"     #%d @ %4d,%4d,%8.1f\n", i, stamp->xss[i], stamp->yss[i], peaks[qs[pcnt-j-1]]); 
            stamp->nss++;
        }
    }
    
    free(qs);
    free(xloc);
    free(yloc);
    free(peaks);
    return 0;
}


/**
 * @brief Compute the mean normalised chi-squared (signal-to-noise-squared) over
 *        a full-frame image using a supplied noise map.
 *
 * @details Iterates over all pixels in the rPixX × rPixY image, skipping pixels
 * that fail the mask selection or have near-zero data values.  For accepted
 * pixels accumulates data[i]² / noise[i]² and divides by the count.  The
 * result in @p nnorm is a dimensionless noise metric; values near 1 indicate
 * a well-normalised difference image.
 *
 * Mask logic (umask / smask pairs select pixel quality):
 *  - good pixels: umask=0,      smask=0xffff
 *  - ok pixels:   umask=0xff,   smask=0x8000
 *  - bad pixels:  umask=0x8000, smask=0
 *
 * @param data     Full-frame data image (e.g. difference image).
 * @param noise    Full-frame per-pixel noise (sigma) image.
 * @param nnorm    Output: mean chi-squared per pixel, or MAXVAL if n < 2.
 * @param nncount  Output: number of valid pixels used in the sum.
 * @param umask    Pixels with this bit set in mRData are excluded.
 * @param smask    Pixels without this bit set in mRData are excluded (0 disables).
 */
void getNoiseStats3(float *data, float *noise, double *nnorm, int *nncount, int umask, int smask) {

    double noiseWeightedSum=0;
    int pixelIdx, validPixelCount=0, maskData;
    float dataValue=0, noiseInv=0;


    for (pixelIdx = rPixX*rPixY; pixelIdx--; ) {
        dataValue = data[pixelIdx];
        maskData = mRData[pixelIdx];

        /*
          fprintf(stderr, "CAW %d %f %f %f %d %d %d %d\n", validPixelCount, noiseWeightedSum, dataValue, noiseInv, maskData,
          ((umask > 0) && (!(maskData & umask))),
          ((smask > 0) &&   (maskData & smask)),
          (dataValue == fillVal) || (fabs(dataValue) <= ZEROVAL));
        */

        if (((umask > 0) && (!(maskData & umask))) ||
            ((smask > 0) &&   (maskData & smask))  ||
            (fabs(dataValue) <= ZEROVAL))
            continue;

        noiseInv  = 1. / noise[pixelIdx];

        validPixelCount    += 1;
        noiseWeightedSum += dataValue*dataValue * noiseInv*noiseInv;
    }
    if (validPixelCount > 1) {
        *nncount = validPixelCount;
        *nnorm   = noiseWeightedSum / (float)validPixelCount;
    }
    else {
        *nncount = validPixelCount;
        *nnorm   = MAXVAL;
    }
    return;
}

/**
 * @brief Compute robust pixel-distribution statistics (mean, median, mode, sigma,
 *        FWHM) for a sub-region of an image using iterative sigma-clipping and
 *        histogram analysis.
 *
 * @details Algorithm (multi-step histogram refinement):
 *
 *   1. **Sampling & bin-width estimation:** Draw HISTOGRAM_SAMPLE_SIZE (100)
 *      random pixels; compute percentiles at HISTOGRAM_LOWER_FRAC (50%) and
 *      HISTOGRAM_UPPER_FRAC (90%) to estimate the dynamic range and bin width.
 *
 *   2. **Sigma-clipping setup:** Call sigma_clip() on all valid pixels to
 *      obtain mean and standard deviation, excluding outliers beyond maxiter
 *      iterations of the clipping loop.
 *
 *   3. **Histogram construction & refinement loop:** Build a 256-bin histogram
 *      and test if bin boundaries fit the data:
 *      - If bins too wide (lower < 1, upper > 255): double bin width, retry.
 *      - If bins too narrow (peak span < 40 bins): halve bin width, retry.
 *      - If bins acceptable: exit loop.
 *      Maximum MAX_SIGMA_CLIP_RETRIES (5) attempts to prevent divergence.
 *
 *   4. **Mode finding:** Locate the histogram peak as the narrowest region
 *      containing ~HISTOGRAM_PEAK_WIDTH_PCT (10%) of pixels; interpolate
 *      bin boundaries using weighted centroid. Store peak location in *mode.
 *
 *   5. **FWHM estimation:** Find 25th and 75th percentiles (i.e., 50% of
 *      pixels around the mode) via cumulative histogram; estimate FWHM as
 *      the difference between these percentiles (a robust surrogate for
 *      the Gaussian sigma of the noise distribution).
 *
 *   Mask filtering uses two bitmask parameters:
 *   - umask (exclude bits): pixels with (mRData[i] & umask) != 0 are skipped.
 *   - smask (select bits): pixels without (mRData[i] & smask) set are skipped.
 *   Good pixels: umask=0x0, smask=0xffff. Semi-bad: umask=0xff, smask=0x8000.
 *
 *   Reference: Alard & Lupton (1998), Appendix; histogram-based estimators
 *   are robust to outliers and provide unbiased FWHM estimates for PSF-convolved
 *   noise (which is non-Gaussian due to correlation).
 *
 * @param data    Flat pixel array of the stamp region (size nPixX × nPixY).
 * @param x0Reg   X-coordinate offset of the stamp within the full region (for mask lookup).
 * @param y0Reg   Y-coordinate offset of the stamp within the full region.
 * @param nPixX   Width of the stamp in pixels.
 * @param nPixY   Height of the stamp in pixels.
 * @param sum     Output: sum of all valid pixel values.
 * @param mean    Output: mean of sigma-clipped pixels.
 * @param median  Output: median (50th percentile) from histogram.
 * @param mode    Output: mode (peak bin) from histogram, interpolated via
 *                weighted centroid of bins near the peak.
 * @param sd      Output: standard deviation (sigma) of sigma-clipped pixels.
 * @param fwhm    Output: full-width half-maximum FWHM estimate (in pixel
 *                value units, not angular FWHM); derived from 75th–25th
 *                percentile range, robust to outliers.
 * @param lfwhm   Output: lower-quartile FWHM estimate (25th percentile).
 * @param umask   Bits to exclude: pixels with (mRData[pixel] & umask) != 0 are skipped.
 * @param smask   Bits to select: pixels without (mRData[pixel] & smask) set are skipped.
 * @param maxiter Maximum number of sigma-clipping iterations in step 2.
 *
 * @return 0 on success; error codes:
 *   - 1: memory allocation failure or max retries exceeded
 *   - 2: insufficient valid pixels (< HISTOGRAM_SAMPLE_SIZE) after masking
 *   - 3: degenerate bin width (zero variance in data)
 *   - 4: too few input pixels (nPixX*nPixY < HISTOGRAM_SAMPLE_SIZE)
 *   - 5: sigma_clip() failed to converge
 */
int getStampStats3(float *data,
                   int x0Reg, int y0Reg, int nPixX, int nPixY,
                   double *sum, double *mean, double *median,
                   double *mode, double *sd, double *fwhm, double *lfwhm,
                   int umask, int smask, int maxiter) {
    /*
      good pixels : umask=0,      smask=0xffff
      ok   pixels : umask=0xff,   smask=0x8000
      bad  pixels : umask=0x8000, smask=0
      
      this version 3 does sigma clipping 
    */
    
    /* this came primarily from Gary Bernstein */
    
    extern int flcomp();
    
    double   binLowerBound,binWidth,maxDensity,modeNormalized;
    double   binSum,weightedBinSum,inverseStdev;
    double   lowerPercentile,upperPercentile,modeBin, pixelValue;
    int      maskData,totalPixelCount;
    int      histogram[256],pixelIdx,rowIdx,randomX,randomY;
    int      binIndex,maxDensityBinIdx=0,refinementAttempts,lowerBinIdx,upperBinIdx, continueRefinement;
    double   absolutePixelSum;
    int      validPixelCount;
    int      randomSeed;
    float    *sdat;
    double   *work;
    
    int      nstat     = HISTOGRAM_SAMPLE_SIZE;
    float    ufstat    = HISTOGRAM_UPPER_FRAC;
    float    mfstat    = HISTOGRAM_LOWER_FRAC;
    
    
    totalPixelCount = nPixX * nPixY;
    if (totalPixelCount < nstat)
        return 4;

    /* fprintf(stderr, "DOINK!  %d %d %f\n", x0Reg, y0Reg, data[0]); */
    if ( !(sdat = (float *)calloc(totalPixelCount, sizeof(float))) ||
         !(work = (double *)calloc(nstat, sizeof(double)))) {
        return (1);
    }
    
    randomSeed   = RNG_SEED_MAGIC;  /* Numerical Recipes convention; triggers re-seeding */
    refinementAttempts  = 0;     /* attempts at the histogram */

    /* ====================================================================
       SAMPLING PHASE: Estimate bin width from percentile range
       ===================================================================== */
    validPixelCount = 0;
    /* pull HISTOGRAM_SAMPLE_SIZE random enough values to estimate required bin sizes */
    /* only do as many calls as there are points in the section! */
    /* NOTE : ignore anything with fillVal or zero */
    for (pixelIdx = 0; (pixelIdx < nstat) && (validPixelCount < totalPixelCount); pixelIdx++, validPixelCount++) {
        randomX = (int)floor(ran1(&randomSeed)*nPixX);
        randomY = (int)floor(ran1(&randomSeed)*nPixY);

        /* here data is size rPixX, rPixY */
        pixelValue = data[randomX+randomY*nPixX];
        /* region is rPixX, rPixY */
        maskData = mRData[(randomX+x0Reg)+(randomY+y0Reg)*rPixX];

        if (((umask > 0) && (!(maskData & umask))) ||
            ((smask > 0) &&   (maskData & smask))  ||
            (fabs(pixelValue) <= ZEROVAL))
            pixelIdx--;
        else
            work[pixelIdx] = pixelValue;
    }
    qsort(work, pixelIdx, sizeof(double), flcomp);
    totalPixelCount = pixelIdx;

    binWidth = (work[(int)(ufstat*totalPixelCount)] - work[(int)(mfstat*totalPixelCount)]) / (float)nstat; /* good estimate */
    binLowerBound    = work[(int)(mfstat*totalPixelCount)] - 128. * binWidth;    /* we use 256 bins */
    
    /*** DO ONLY ONCE! ***/
    validPixelCount = 0;
    for (rowIdx = 0; rowIdx < nPixY; rowIdx++) {
        for (pixelIdx = 0; pixelIdx < nPixX; pixelIdx++) {
            pixelValue = data[pixelIdx+rowIdx*nPixX];
            maskData = mRData[(pixelIdx+x0Reg)+(rowIdx+y0Reg)*rPixX];

            /* fprintf(stderr, "BIGT %d %d %f : %d %d %d\n", pixelIdx, rowIdx, pixelValue, pixelIdx+x0Reg, rowIdx+y0Reg, maskData); */
            /* looks like it works.  that is, the image pixel values
               pixelIdx+x0Reg, rowIdx+y0Reg are the same as printed pixel values
               pixelIdx,rowIdx,pixelValue.  mask should work similarly */

            if (((umask > 0) && (!(maskData & umask))) ||
                ((smask > 0) &&   (maskData & smask))  ||
                (fabs(pixelValue) <= ZEROVAL))
                continue;

            if (pixelValue*0.0 != 0.0) {
                mRData[(pixelIdx+x0Reg)+(rowIdx+y0Reg)*rPixX] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                /* fprintf(stderr, "OUCH %d %d %f %d\n", pixelIdx+x0Reg, rowIdx+y0Reg, pixelValue, maskData); */
                continue;
            }

            /* good pixels pass, so sigma clipping part here */
            sdat[validPixelCount++] = pixelValue;
        }
    }

    /* ====================================================================
       SIGMA-CLIPPING PHASE: Obtain mean & std dev before histogram
       ===================================================================== */
    if (sigma_clip(sdat, validPixelCount, mean, sd, maxiter)) {
        free(sdat);
        free(work);
        return 5;
    }
    free(sdat);

    /* save some speed */
    inverseStdev = 1. / (*sd);

    /* ====================================================================
       HISTOGRAM REFINEMENT LOOP: Adjust bin width until bounds are acceptable
       ===================================================================== */
    continueRefinement = 1;
    while (continueRefinement) {

        if (refinementAttempts >= MAX_SIGMA_CLIP_RETRIES) {
            /* too many attempts here - print message and exit*/
            /* DD fprintf(stderr, "     WARNING: 5 failed iterations in getStampStats2\n"); */
            free(work);
            return 1;
        }

        for (pixelIdx=0; pixelIdx<HISTOGRAM_NUM_BINS; histogram[pixelIdx++]=0);
        
        /* rezero sums if repeating */
        absolutePixelSum = binSum = weightedBinSum = 0.;
        validPixelCount = 0;

        for (rowIdx = 0; rowIdx < nPixY; rowIdx++) {
            for (pixelIdx = 0; pixelIdx < nPixX; pixelIdx++) {
                pixelValue = data[pixelIdx+rowIdx*nPixX];
                maskData = mRData[(pixelIdx+x0Reg)+(rowIdx+y0Reg)*rPixX];

                if (((umask > 0) && (!(maskData & umask))) ||
                    ((smask > 0) &&   (maskData & smask))  ||
                    (fabs(pixelValue) <= ZEROVAL))
                    continue;

                if (pixelValue*0.0 != 0.0) {
                    mRData[(pixelIdx+x0Reg)+(rowIdx+y0Reg)*rPixX] |= (FLAG_INPUT_ISBAD | FLAG_ISNAN);
                    continue;
                }

                /* final sigma cut */
                /* reject both high and low here */
                if ((fabs(pixelValue - (*mean)) * inverseStdev) > statSig)
                    continue;

                binIndex = floor( (pixelValue-binLowerBound)/binWidth ) + 1;
                binIndex = (binIndex < 0 ? 0 : binIndex);
                binIndex = (binIndex > 255 ? 255 : binIndex);

                /* NOTE : the zero index is way overweighted since it contains everything
                   below binLowerBound.
                */
                histogram[binIndex]++;

                /* absolutePixelSum is the sum of absolute value of pixels in the image */
                absolutePixelSum += fabs(pixelValue);
                /*
                  This is the number of pixels we are working with, not totalPixelCount which
                  includes any zeroed or masked pixels.
                */
                validPixelCount += 1;

            }
        }

        /* get out of dodge! */
        if (validPixelCount == 0) {
            /* DD fprintf(stderr, "     WARNING: no stampStats2 data\n"); */
            *mode = *median = work[(int)(mfstat*totalPixelCount)];
            *fwhm = *lfwhm = 0.;
            free(work);
            return 2;
        }

        /* quit here if the bins are degenerate */
        if (binWidth == 0.) {
            /* DD fprintf(stderr, "     WARNING: no variation in stampStats2 data\n"); */
            *mode = *median = work[(int)(mfstat*totalPixelCount)];
            *fwhm = *lfwhm = 0.;
            free(work);
            return 3;
        }

        /* ====================================================================
           MODE FINDING: Locate peak as narrowest region containing ~10% of data
           ===================================================================== */
        /* find the mode - find narrowest region which holds ~10% of points*/
        binSum = maxDensity = 0.;
        for (lowerBinIdx = upperBinIdx = 1; upperBinIdx < 255; binSum -= histogram[lowerBinIdx++]) {
            while ( (binSum < validPixelCount * HISTOGRAM_PEAK_WIDTH_PCT) && (upperBinIdx < 255) )
                binSum += histogram[upperBinIdx++];

            if (binSum / (upperBinIdx-lowerBinIdx) > maxDensity) {
                maxDensity = binSum / (upperBinIdx-lowerBinIdx);
                maxDensityBinIdx = lowerBinIdx;
            }
        }
        /* if it never got assigned... */
        if (maxDensityBinIdx < 0 || maxDensityBinIdx > 255)
            maxDensityBinIdx = 0;


        /* try to interpolate between bins somewhat by finding weighted
           mean of the bins in the peak*/
        /* NOTE : need <= here and above in case validPixelCount is small (< 10) */
        weightedBinSum = binSum = 0.;
        for (pixelIdx = maxDensityBinIdx; (binSum < validPixelCount/10.) && (pixelIdx < 255); pixelIdx++) {
            binSum  += histogram[pixelIdx];
            weightedBinSum += pixelIdx*histogram[pixelIdx];
        }
        modeBin = weightedBinSum / binSum + 0.5; /*add 0.5 to give middle of bin*/
        *mode = binLowerBound + binWidth * (modeBin - 1.);

        /*find the percentile of the mode*/
        maxDensityBinIdx = floor(modeBin);
        for (pixelIdx = 0, binSum = 0.; pixelIdx < maxDensityBinIdx; binSum += histogram[pixelIdx++]);
        binSum += histogram[maxDensityBinIdx] * (modeBin-maxDensityBinIdx); /*interpolate fractional bin*/
        binSum /= validPixelCount;
        modeNormalized=binSum;

        /* ====================================================================
           FWHM ESTIMATION: Find 25th–75th percentile range (robust noise σ estimate)
           ===================================================================== */
        /* find the region around mode containing half the "noise" points,
           e.g. assume that the mode is 50th percentile of the noise */
        lowerPercentile = validPixelCount * (0.5 - HISTOGRAM_NOISE_HALF_PCT);  /* 25th percentile */
        upperPercentile = validPixelCount * (0.5 + HISTOGRAM_NOISE_HALF_PCT);  /* 75th percentile */
        binSum = 0.;
        for (pixelIdx = 0; binSum < lowerPercentile; binSum += histogram[pixelIdx++]);
        lowerPercentile = pixelIdx - (binSum - lowerPercentile) / histogram[pixelIdx-1];
        for ( ; binSum < upperPercentile; binSum += histogram[pixelIdx++]);
        upperPercentile = pixelIdx - (binSum - upperPercentile) / histogram[pixelIdx-1];

        /*now a few checks to make sure the histogram bins were chosen well*/
        if ( (lowerPercentile < 1.) || (upperPercentile > 255.) ) {
            /*make the bins wider, about same center*/
            binLowerBound -= 128. * binWidth;
            binWidth *= 2.;
            refinementAttempts++;
            continueRefinement = 1;
        } else if ( (upperPercentile-lowerPercentile) < 40. ) {
            /*make bins narrower for better precision*/
            binWidth /= 3.;
            binLowerBound = *mode - 128. * binWidth;
            refinementAttempts++;
            continueRefinement=1;;
        } else
            continueRefinement = 0;
        
        /* end of re-histogramming loop */
    }
    
    *sum = absolutePixelSum;

    /* calculate the noise sd based on this distribution width
       the numerical constant converts the width to sd based on a
       gaussian noise distribution
    */
    *fwhm = binWidth * (upperPercentile - lowerPercentile) / 1.35;

    /*find the median*/
    for (pixelIdx = 0, binSum = 0; binSum < validPixelCount/2.; binSum += histogram[pixelIdx++]);
    *median = pixelIdx - (binSum - validPixelCount/2.) / histogram[pixelIdx-1];

    /* lower-quartile result - scale to sigma: */
    *lfwhm = binWidth * (*median - lowerPercentile) * 2. / 1.35;

    *median = binLowerBound + binWidth*(*median-1.);

    free(work);
    return 0;
}

/**
 * @brief Iteratively sigma-clip a float array to compute a robust mean and
 *        standard deviation.
 *
 * @details Repeatedly computes the mean and standard deviation of all
 * unmasked points, then masks any point more than statSig standard deviations
 * from the current mean.  Convergence is declared when no additional points are
 * rejected or @p maxiter iterations have been completed.  Passing maxiter=0
 * computes the unclipped mean and standard deviation.
 *
 * @param data     Input float array of length @p count.
 * @param count    Number of elements in @p data.
 * @param mean     Output: sigma-clipped mean.
 * @param stdev    Output: sigma-clipped standard deviation.
 * @param maxiter  Maximum number of rejection iterations (0 = no clipping).
 * @return 0 on success; 1 if count == 0; 2 if all points rejected; 3 if only
 *         one point remains.
 */
int sigma_clip(float *data, int count, double *mean, double *stdev, int maxiter) {
    int acceptedCount, newAcceptedCount, dataIdx, iterationCount;
    char *clipMask;
    double inverseStdev;
    float dataValue;

    if (count == 0) {
        *mean  = 0;
        *stdev = MAXVAL;
        return 1;
    }

    clipMask = (char *)calloc(count, sizeof(char));
    /*for (dataIdx=0; dataIdx<count; dataIdx++) clipMask[dataIdx] = 0;*/

    acceptedCount = 0;
    newAcceptedCount = count;
    iterationCount = 0;

    while ((newAcceptedCount != acceptedCount) && (iterationCount < maxiter)) {
        acceptedCount = newAcceptedCount;

        *mean  = 0;
        *stdev = 0;
        for (dataIdx=0; dataIdx<count; dataIdx++) {
            if (!(clipMask[dataIdx])) {
                dataValue       = data[dataIdx];
                *mean  += dataValue;
                *stdev += dataValue*dataValue;
            }
        }

        if (newAcceptedCount > 0)
            *mean /= newAcceptedCount;
        else {
            *mean  = 0;
            *stdev = MAXVAL;
            free(clipMask);
            return 2;
        }

        if (newAcceptedCount > 1) {
            *stdev = *stdev - newAcceptedCount * (*mean) * (*mean);
            *stdev = sqrt(*stdev / (double)(newAcceptedCount - 1));
        }
        else {
            *stdev = MAXVAL;
            free(clipMask);
            return 3;
        }

        newAcceptedCount   = 0;
        inverseStdev = 1. / (*stdev);
        for (dataIdx=0; dataIdx<count; dataIdx++) {
            if (!(clipMask[dataIdx])) {
                /* reject high and low outliers */
                if ((fabs(data[dataIdx] - (*mean)) * inverseStdev) > statSig) {
                    clipMask[dataIdx] = 1;
                }
                else {
                    newAcceptedCount++;
                }
            }
        }
        iterationCount += 1;
        /*fprintf(stderr, "%d %d %f %f\n", acceptedCount, newAcceptedCount, (*mean), (*stdev));*/
    }
    /*fprintf(stderr, "\n");*/
    free(clipMask);
    return 0;
}

/**
 * @brief Compute a spatially-smoothed noise or mean image using a sliding box
 *        window.
 *
 * @details For each pixel (i, j), collects the values of all unmasked
 * neighbours within a (2*size+1)×(2*size+1) box, calls sigma_clip() with
 * maxiter=0 (no clipping), and stores either the box mean (doavg=1, useful for
 * a smoothed noise map) or the box standard deviation (doavg=0, useful for an
 * empirical RMS noise image from the difference image).  Boundary pixels are
 * handled by simply skipping out-of-bounds neighbours.
 *
 * @param image    Input nx×ny float image.
 * @param mask     Integer mask array (same size); pixels with (mask[k] & maskval)
 *                 non-zero are excluded from the box statistics.
 * @param nx       Image width.
 * @param ny       Image height.
 * @param size     Half-width of the sliding box (box is (2*size+1) pixels wide).
 * @param maskval  Mask bitmask: pixels with this bit set are excluded.
 * @param doavg    1 to return the box mean; 0 to return the box standard deviation.
 * @return Pointer to a newly allocated nx×ny output float array, or NULL on
 *         allocation failure.
 */
float *calculateAvgNoise(float *image, int *mask, int nx, int ny, int size, int maskval, int doavg) {
    int pixelX, pixelY, neighborColIdx, neighborRowIdx, validPixelCount;
    float *data, *outim;
    double mean, stdev;


    if ( !(data = (float *) calloc(size * size, sizeof(float))))
        return (NULL);
    if ( !(outim = (float *) calloc(nx * ny, sizeof(float))))
        return (NULL);

    for (pixelY = ny; pixelY--; ) {
        for (pixelX = nx; pixelX--; ) {

            /* take your average! */
            validPixelCount = 0;
            for (neighborRowIdx = pixelY - size; neighborRowIdx <= pixelY + size; neighborRowIdx++) {
                if ((neighborRowIdx < 0) || neighborRowIdx >= ny)
                    continue;

                for (neighborColIdx = pixelX - size; neighborColIdx <= pixelX + size; neighborColIdx++) {
                    if ((neighborColIdx < 0) || neighborColIdx >= nx)
                        continue;

                    if (!(mask[neighborColIdx+neighborRowIdx*nx] & maskval))
                        data[validPixelCount++] = image[neighborColIdx+neighborRowIdx*nx];
                }
            }
            /* we'll do no clipping */
            sigma_clip(data, validPixelCount, &mean, &stdev, 0);
            outim[pixelX+pixelY*nx] = doavg ? mean : stdev;
        }
    }
    free(data);
    return outim;
}


/**
 * @brief Release all dynamically allocated sub-arrays within an array of
 *        stamp_struct objects.
 *
 * @details Frees all dynamically allocated arrays within each stamp_struct:
 * vectors, mat, krefArea, scprod, xss, and yss.
 *
 * The new contiguous allocation strategy (alloc_vector_array, alloc_matrix_contiguous)
 * means vectors and mat each require only 2 free() calls total instead of N+1.
 *
 * @param stamps   Array of stamp_struct objects whose sub-arrays are to be freed.
 * @param nStamps  Number of stamps in the array.
 */
void freeStampMem(stamp_struct *stamps, int nStamps) {
    int stampIdx;

    if (stamps) {
        for (stampIdx = 0; stampIdx < nStamps; stampIdx++) {
            /* Free contiguous vector array (2 frees: pointer array + data block) */
            if (stamps[stampIdx].vectors) free_vector_array(stamps[stampIdx].vectors);

            /* Free contiguous matrix (2 frees: pointer array + data block) */
            if (stamps[stampIdx].mat) free_matrix_contiguous(stamps[stampIdx].mat);

            /* Free simple arrays */
            if (stamps[stampIdx].krefArea) free(stamps[stampIdx].krefArea);
            if (stamps[stampIdx].scprod) free(stamps[stampIdx].scprod);
            if (stamps[stampIdx].xss) free(stamps[stampIdx].xss);
            if (stamps[stampIdx].yss) free(stamps[stampIdx].yss);
        }
    }
}

/**
 * @brief Construct a per-pixel noise (sigma) image from an input image using a
 *        Poisson-plus-readout noise model.
 *
 * @details For each pixel computes:
 *   noise² = |iData[i]| · invGain + quad²
 * where invGain = 1/gain converts detected electrons to ADU, and quad is the
 * quadrature readout noise term in ADU.  Returns a noise image (not variance)
 * suitable for chi-squared weighting.
 *
 * @param iData    Full-frame input image in ADU.
 * @param invGain  Reciprocal of the CCD gain (electrons per ADU)^{-1}.
 * @param quad     Readout noise in ADU (added in quadrature).
 * @return Pointer to a newly allocated rPixX×rPixY noise image, or NULL on
 *         allocation failure.
 */
float *makeNoiseImage4(float *iData, float invGain, float quad) {
    
    int    i;
    double qquad;
    float  *nData=NULL;
    
    if ( !(nData = (float *)calloc(rPixX*rPixY, sizeof(float))))
        return NULL;
    
    qquad = quad * quad;
    
    for (i = rPixX*rPixY; i--; ) 
        nData[i] = fabs(iData[i])*invGain + qquad;
    
    return nData;
}

/**
 * @brief Read kernel configuration parameters from an existing kernel FITS
 *        file header, overriding defaults and command-line options.
 *
 * @details Opens the kernel FITS file, verifies the KERINFO keyword, then reads
 * NGAUSS, FWKERN, CKORDER, BGORDER, PHOTNORM, and per-Gaussian DGAUSSn /
 * SGAUSSn keywords from the final binary table HDU.  The sigma values read as
 * FWHM-like quantities are converted to the internal representation
 * sigma_gauss[i] = 1/(2·σ²).  This allows a pre-computed kernel to be applied
 * to a new pair of images without repeating the fit.
 *
 * @param kimage  Path to the kernel FITS file produced by a prior HOTPANTS run.
 */
void getKernelInfo(char *kimage) {
    
    fitsfile *kPtr;
    int i, existsTable, status = 0;
    char hKeyword[1024];
    
    /* open the input kernel image */
    if ( fits_open_file(&kPtr, kimage, 0, &status) )
        printError(status);
    
    /* required keyword in primary HDU */
    if ( fits_read_key_log(kPtr, "KERINFO", &existsTable, NULL, &status) )
        printError(status);
    
    if (!(existsTable)) {
        fits_close_file(kPtr, &status);
        fprintf(stderr, "This image does not appear to contain a kernel table, exiting...\n");
        exit(1);
    }
    
    /* move to binary kernel table... */
    if ( fits_get_num_hdus(kPtr, &existsTable, &status) ||
         fits_movabs_hdu(kPtr, existsTable, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "NGAUSS", &ngauss, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "FWKERN", &fwKernel, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "CKORDER", &kerOrder, NULL, &status) ||
         fits_read_key(kPtr, TINT,    "BGORDER", &bgOrder, NULL, &status) )
        printError(status);
    
    deg_fixe    = (int *)realloc(deg_fixe,      ngauss*sizeof(int));
    sigma_gauss = (float *)realloc(sigma_gauss, ngauss*sizeof(float));
    
    /* this took a while to figure out! */
    photNormalize = (char *)malloc(1*sizeof(char));
    
    snprintf(hKeyword, sizeof(hKeyword), "PHOTNORM");
    if (fits_read_key(kPtr, TSTRING, hKeyword, photNormalize, NULL, &status))
        printError(status);
    
    /* read kernel gaussian info */
    for (i = 0; i < ngauss; i++) {
        snprintf(hKeyword, sizeof(hKeyword), "DGAUSS%d", i+1);
        if (fits_read_key(kPtr, TINT, hKeyword, &deg_fixe[i], NULL, &status))
            printError(status);
        snprintf(hKeyword, sizeof(hKeyword), "SGAUSS%d", i+1);
        if (fits_read_key(kPtr, TFLOAT, hKeyword, &sigma_gauss[i], NULL, &status))
            printError(status);
        
        /* important! */
        sigma_gauss[i] = (1.0/(2.0*sigma_gauss[i]*sigma_gauss[i]));
    }
    
    if (fits_close_file(kPtr, &status) )
        printError(status);
    
    return;
}


/**
 * @brief Read the kernel solution for one image region from a pre-computed
 *        kernel FITS file.
 *
 * @details Reads from the FITS header the spatial extent of the region
 * (REGIONnn keyword, FITS 1-indexed, converted to 0-indexed internally), the
 * convolution direction (CONVOLnn), quality-control statistics (SSSIGnn,
 * SSSSCATnn, FSIGnn, FSCATnn, NSCALOnn), and the kernel solution coefficients
 * from the binary table column nRegion+1.  Depending on the convolution
 * direction, the coefficients are stored in *tKerSol (template convolved) or
 * *iKerSol (image convolved).
 *
 * @param kimage               Path to the kernel FITS file.
 * @param nRegion              0-based region index.
 * @param tKerSol              In/out: pointer to the template kernel solution
 *                             vector; reallocated to nCompTotal+1 doubles if
 *                             the template was convolved.
 * @param iKerSol              In/out: pointer to the image kernel solution
 *                             vector; reallocated similarly.
 * @param rXMin                Output: left pixel bound of the region (0-indexed).
 * @param rXMax                Output: right pixel bound.
 * @param rYMin                Output: bottom pixel bound.
 * @param rYMax                Output: top pixel bound.
 * @param meansigSubstamps     Output: mean substamp residual sigma from the fit.
 * @param scatterSubstamps     Output: scatter of substamp residuals from the fit.
 * @param meansigSubstampsF    Output: mean substamp residual sigma from the
 *                             final iteration.
 * @param scatterSubstampsF    Output: scatter from the final iteration.
 * @param diffrat              Output: photometric scale factor ratio (NSCALO);
 *                             defaults to 1 if the keyword is absent.
 * @param NskippedSubstamps    Output: number of skipped substamps (unused here,
 *                             kept for API consistency).
 */
void readKernel(char *kimage, int nRegion, double **tKerSol, double **iKerSol,
                int *rXMin, int *rXMax, int *rYMin, int *rYMax,
                double *meansigSubstamps, double *scatterSubstamps,
                double *meansigSubstampsF, double *scatterSubstampsF,
                double *diffrat, int *NskippedSubstamps) {
    
    fitsfile *kPtr;
    int status = 0;
    char hKeyword[1024], hInfo[1024];
    
    /* open the input kernel image */
    if ( fits_open_file(&kPtr, kimage, 0, &status) )
        printError(status);
    
    /* grab stuff for this region */
    snprintf(hKeyword, sizeof(hKeyword), "REGION%02d", nRegion);
    if (fits_read_key(kPtr, TSTRING, hKeyword, &hInfo, NULL, &status))
        printError(status);
    
    /* get extent of region */
    if (sscanf(hInfo, "[%d:%d,%d:%d]", rXMin, rXMax, rYMin, rYMax) != 4) {
        fprintf(stderr, "Problem with region %d (%s), exiting...\n", nRegion, hInfo);
        exit(1);
    }
    /* fits indexing starts at 1, code at 0 */
    *rXMin -= 1;
    *rXMax -= 1;
    *rYMin -= 1;
    *rYMax -= 1;
    
    /* which way to convolve */
    snprintf(hKeyword, sizeof(hKeyword), "CONVOL%02d", nRegion);
    if (fits_read_key(kPtr, TSTRING, hKeyword, &hInfo, NULL, &status))
        printError(status);
    
    /* copy quality control stuff: mean sigma, scatter, # substamps skipped */
    snprintf(hKeyword, sizeof(hKeyword), "SSSIG%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, meansigSubstamps, NULL, &status))
        printError(status);
    
    snprintf(hKeyword, sizeof(hKeyword), "SSSCAT%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, scatterSubstamps, NULL, &status))
        printError(status);
    
    snprintf(hKeyword, sizeof(hKeyword), "FSIG%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, meansigSubstampsF, NULL, &status))
        printError(status);
    
    snprintf(hKeyword, sizeof(hKeyword), "FSCAT%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, scatterSubstampsF, NULL, &status))
        printError(status);
    
    /* sometimes does not exist */
    snprintf(hKeyword, sizeof(hKeyword), "NSCALO%02d", nRegion);
    if (fits_read_key(kPtr, TDOUBLE, hKeyword, diffrat, NULL, &status)) {
        *diffrat = 1;
        status = 0;
    }
    
    if (strncmp(hInfo, "TEMPLATE", 8)==0) {
        forceConvolve = "t";
        *tKerSol = (double *)realloc(*tKerSol, (nCompTotal+1)*sizeof(double));
        
        fits_get_kernel_btbl(kPtr, &(*tKerSol), nRegion);
    }
    else if (strncmp(hInfo, "IMAGE", 5)==0) {
        forceConvolve = "i";
        *iKerSol = (double *)realloc(*iKerSol, (nCompTotal+1)*sizeof(double));
        
        fits_get_kernel_btbl(kPtr, &(*iKerSol), nRegion);
    }
    
    if (fits_close_file(kPtr, &status) )
        printError(status);
    
    return;
}

/**
 * @brief Read one region's kernel solution coefficients from the binary table
 *        extension of a kernel FITS file.
 *
 * @details Moves to the last HDU of the open FITS file (the binary table
 * containing all region solutions), then reads nCompTotal+1 double-precision
 * values from column nRegion+1 into *kernelSol, zeroing the array first.
 * The table layout stores one column per image region and one row per
 * kernel coefficient.
 *
 * @param kPtr       Open cfitsio fitsfile pointer (read mode).
 * @param kernelSol  In/out: pointer to a pre-allocated (nCompTotal+1) double
 *                   array; filled with the kernel solution on return.
 * @param nRegion    0-based region index; reads from table column nRegion+1.
 */
void fits_get_kernel_btbl(fitsfile *kPtr, double **kernelSol, int nRegion) {
    int status=0, existsTable;
    
    /* move to binary kernel table... */
    if ( fits_get_num_hdus(kPtr, &existsTable, &status) ||
         fits_movabs_hdu(kPtr, existsTable, NULL, &status) )
        printError(status);
    
    memset(*kernelSol, 0, (nCompTotal+1)*sizeof(double));
    if (fits_read_col(kPtr, TDOUBLE, nRegion+1, 1, 1, (nCompTotal+1), 0, *kernelSol, 0, &status))
        printError(status);
    
    /*
      fprintf(stderr, "OK %d, %f %f %f %f %f %f\n", nRegion, *kernelSol[0], *kernelSol[1],
      *kernelSol[192-1], *kernelSol[293-1],
      *kernelSol[294-1], *kernelSol[298-1]);
      */
    
    return;
}


/**
 * @brief Dilate the FLAG_INPUT_ISBAD mask by @p width pixels in all directions.
 *
 * @details For every pixel already flagged FLAG_INPUT_ISBAD, sets FLAG_OK_CONV
 * on all neighbours within a (width/2) × (width/2) box (if those neighbours
 * are not themselves flagged as bad).  This marks pixels adjacent to bad pixels
 * as "convolved OK but near bad data", so that downstream processing can apply
 * appropriate caution.  Called after makeInputMask() to propagate the mask by
 * hwKernel × kfSpreadMask1 pixels, ensuring that convolution artefacts near
 * bad pixels are correctly identified in the output.
 *
 * @param mData  Full-frame integer mask array (rPixX × rPixY); modified in place.
 * @param width  Dilation full-width in pixels; no-op if <= 0.
 */
void spreadMask(int *mData, int width) {

    int imageCol, imageRow, horzOffset, vertOffset, ii, jj, w2;

    /* nothing to do! */
    if (width <= 0) {
        return;
    }
    w2 = width/2;
    for (imageRow = 0; imageRow < rPixY; imageRow++) {
        for (imageCol = 0; imageCol < rPixX; imageCol++) {
            if (mData[imageCol+rPixX*imageRow] & FLAG_INPUT_ISBAD) {
                for (horzOffset = -w2; horzOffset <= w2; horzOffset++) {
                    ii = imageCol + horzOffset;
                    if (ii < 0 || ii >= rPixX)
                        continue;

                    for (vertOffset = -w2; vertOffset <= w2; vertOffset++) {
                        jj = imageRow + vertOffset;
                        if (jj < 0 || jj >= rPixY)
                            continue;

                        mData[ii+rPixX*jj] |= FLAG_OK_CONV * (!(mData[ii+rPixX*jj] & FLAG_INPUT_ISBAD));
                    }
                }
            }
        }
    }
    return;
}

/**
 * @brief Build the combined input bad-pixel mask from template and image data,
 *        then dilate it by the kernel half-width.
 *
 * @details Marks pixels in mData with FLAG_INPUT_ISBAD and the appropriate
 * reason flag for three conditions:
 *  - Pixel equals fillVal in either image (FLAG_BAD_PIXVAL).
 *  - Pixel exceeds the upper threshold tUThresh or iUThresh (FLAG_SAT_PIXEL).
 *  - Pixel falls below the lower threshold tLThresh or iLThresh (FLAG_LOW_PIXEL).
 * After setting the initial mask, calls spreadMask() with a dilation of
 * hwKernel × kfSpreadMask1 pixels so that pixels whose convolution footprint
 * overlaps a bad pixel are correctly flagged in the output difference image.
 *
 * @param tData  Full-frame template image.
 * @param iData  Full-frame science image.
 * @param mData  Full-frame integer mask array; updated in place.
 */
void makeInputMask(float *tData, float *iData, int *mData) {

    int pixelIdx;

    for (pixelIdx = rPixX*rPixY; pixelIdx--; ){
        mData[pixelIdx] |= (FLAG_INPUT_ISBAD | FLAG_BAD_PIXVAL) * (tData[pixelIdx] == fillVal  || iData[pixelIdx] == fillVal);
        mData[pixelIdx] |= (FLAG_INPUT_ISBAD | FLAG_SAT_PIXEL)  * (tData[pixelIdx] >= tUThresh || iData[pixelIdx] >= iUThresh);
        mData[pixelIdx] |= (FLAG_INPUT_ISBAD | FLAG_LOW_PIXEL)  * (tData[pixelIdx] <= tLThresh || iData[pixelIdx] <= iLThresh);
    }
    
    spreadMask(mData, (int)(hwKernel*kfSpreadMask1));
    
    /* mask has value 0 for good pixels in the difference image */
    return;
}

/**
 * @brief Copy non-structural FITS header keywords from one file to another,
 *        excluding NAXIS keywords, EXTEND, and boilerplate COMMENT cards.
 *
 * @details Reads nkeys from the source HDU and writes each card to the output
 * HDU, skipping the first (4 + naxis) structural keywords and any EXTEND or
 * standard FITS boilerplate COMMENT cards.  Used when creating output FITS
 * images that should inherit the WCS and instrument keywords of an input image.
 *
 * @param iPtr    Open cfitsio fitsfile pointer to the source HDU.
 * @param oPtr    Open cfitsio fitsfile pointer to the destination HDU.
 * @param status  cfitsio status code; passed through and returned.
 * @return The cfitsio status code after all operations.
 */
int hp_fits_copy_header(fitsfile *iPtr, fitsfile *oPtr, int *status) {
#define FSTRNCMP(a,b,n)  ((a)[0]<(b)[0]?-1:(a)[0]>(b)[0]?1:strncmp((a),(b),(n)))   
    int nkeys, i;
    long naxis;
    char card[SCRLEN];
    
    if (fits_get_hdrspace(iPtr, &nkeys, NULL, status) ||
        fits_read_key_lng(iPtr, "NAXIS", &naxis, NULL, status) )
        return *status;
    
    /* copy remaining keywords, excluding NAXIS?, EXTEND, and reference COMMENT keywords */
    for (i = 4 + naxis; i <= nkeys; i++) {
        if (fits_read_record(iPtr, i, card, status))
            break;
        
        if (FSTRNCMP(card, "EXTEND  ", 8) &&
            FSTRNCMP(card, "COMMENT   FITS (Flexible Image Transport System) format is", 58) &&
            FSTRNCMP(card, "COMMENT   and Astrophysics', volume 376, page 3", 47) ) {
            if (fits_write_record(oPtr, card, status))
                return *status;
        }
    }
    return *status;
}

/**
 * @brief Clamp float pixel values to the valid range for 16-bit integer FITS
 *        output and flag overflowing pixels.
 *
 * @details If @p makeShort is non-zero, computes the representable range
 * [−32768·bScale+bZero, 32767·bScale+bZero] and sets any out-of-range pixel
 * to the boundary value, simultaneously flagging it FLAG_OUTPUT_ISBAD in
 * mRData.  This guards against photometric rescaling driving pixel values
 * beyond the BITPIX=16 dynamic range.
 *
 * @param data       Float pixel array of length @p npix; modified in place.
 * @param npix       Number of pixels.
 * @param bZero      FITS BZERO keyword value.
 * @param bScale     FITS BSCALE keyword value.
 * @param makeShort  If non-zero, apply short-integer range clamping.
 */
void hp_fits_correct_data(float *data, int npix, float bZero, float bScale, int makeShort) {
    
    float maxval=1e30, minval=-1e30;
    int   i;
    float *dptr;
    int   *mptr;
    
    /*
      BUYER BEWARE : the photometric rescaling of the kernel can take
      an innocent amount of flux and drive it higher than the allowed
      short maximum value.  we need to check for this here
      
      %%% AND WATCH OUT FOR BSCALE MADNESS... %%%
      
    */
    
    dptr = data;
    mptr = mRData;
    
    if (makeShort) {
        maxval =  32767. * bScale + bZero;
        minval = -32768. * bScale + bZero;
        
        for (i = 0; i < npix; i++, dptr++, mptr++) {
            if (*dptr > maxval) {
                *dptr  = maxval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
            else if (*dptr < minval) {
                *dptr  = minval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
        }
    }
    /* extra check for NaN, probably not necessary */
    
    /*
      dptr = data;
      for (i = 0; i < npix; i++, dptr++) {
      if (*dptr*0 != 0)
      *dptr = fillVal;
      }
    */
    
    return;
}

/**
 * @brief Clamp integer pixel values to the valid range for 16-bit integer FITS
 *        output and flag overflowing pixels.
 *
 * @details Integer variant of hp_fits_correct_data() applied to mask or integer
 * output images.  Behaviour and parameters are identical to the float version
 * except that @p data is an int array.
 *
 * @param data       Integer pixel array of length @p npix; modified in place.
 * @param npix       Number of pixels.
 * @param bZero      FITS BZERO keyword value.
 * @param bScale     FITS BSCALE keyword value.
 * @param makeShort  If non-zero, apply short-integer range clamping.
 */
void hp_fits_correct_data_int(int *data, int npix, float bZero, float bScale, int makeShort) {
    
    float maxval=1e30, minval=-1e30;
    int   i;
    int  *dptr, *mptr;
    
    /*
      BUYER BEWARE : the photometric rescaling of the kernel can take
      an innocent amount of flux and drive it higher than the allowed
      short maximum value.  we need to check for this here
      
      %%% AND WATCH OUT FOR BSCALE MADNESS... %%%
      
    */
    
    dptr = data;
    mptr = mRData;
    
    if (makeShort) {
        maxval =  32767. * bScale + bZero;
        minval = -32768. * bScale + bZero;
        
        for (i = 0; i < npix; i++, dptr++) {
            if (*dptr > maxval) {
                *dptr  = maxval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
            else if (*dptr < minval) {
                *dptr  = minval;
                *mptr |= FLAG_OUTPUT_ISBAD;
            }
        }
    }
    return;
}

/**
 * @brief Write a rectangular sub-region of a float array to a FITS image,
 *        clamping values to the 16-bit range if required.
 *
 * @details First calls hp_fits_correct_data() to clamp the full rPixX×rPixY
 * buffer in place, then writes rows fpixelY..lpixelY, columns fpixelX..lpixelX
 * to the FITS file using fits_write_subset_flt(), skipping xArrayLo columns
 * at the beginning of each row.  This allows writing a tiled sub-region of the
 * output image without copying data to a separate buffer.
 *
 * @param fptr      Open cfitsio fitsfile pointer (write mode).
 * @param group     Group parameter (0 for primary image).
 * @param naxis     Number of image axes.
 * @param naxes     Array of axis sizes.
 * @param data      Full rPixX×rPixY float image; values clamped in place.
 * @param status    cfitsio status code; passed through.
 * @param makeShort If non-zero, clamp to 16-bit short range.
 * @param bZero     FITS BZERO.
 * @param bScale    FITS BSCALE.
 * @param fpixelX   First pixel column to write (1-indexed FITS convention).
 * @param fpixelY   First pixel row to write.
 * @param lpixelX   Last pixel column to write.
 * @param lpixelY   Last pixel row to write.
 * @param xArrayLo  Column offset within the @p data row for the first pixel.
 * @param yArrayLo  Row offset within @p data for the first row.
 * @return cfitsio status code.
 */
int hp_fits_write_subset(fitsfile *fptr, long group, long naxis, long *naxes,
                         float *data, int *status, int makeShort,
                         float bZero, float bScale,
                         int fpixelX, int fpixelY, int lpixelX, int lpixelY, int xArrayLo, int yArrayLo) {
    
    int   pixX, pixY, y, x;
    float *dptr;
    long  fpixel[2], lpixel[2];
    
    hp_fits_correct_data(data, (rPixX*rPixY), bZero, bScale, makeShort);
    
    pixX = lpixelX - fpixelX + 1;
    pixY = lpixelY - fpixelY + 1;
    
    fpixel[0] = fpixelX;
    lpixel[0] = lpixelX;
    
    dptr = data;
    /* get to the first row */
    for (y = 0; y < yArrayLo; y++)
        for (x = 0; x < rPixX; x++, dptr++);
    
    for (y = 0; y < pixY; y++) {
        
        fpixel[1] = fpixelY + y;
        lpixel[1] = fpixel[1];
        
        /* get to first pixel in row to write */
        for (x = 0; x < xArrayLo; x++, dptr++);
        
        /* this works! */
        fits_write_subset_flt(fptr, group, naxis, naxes, fpixel, lpixel, dptr, status);
        
        /* clear the row, as fits_write_subset does not increment dptr */
        for (x = xArrayLo; x < rPixX; x++, dptr++);
    }
    return *status;
}

/**
 * @brief Write a rectangular sub-region of an integer array to a FITS image,
 *        clamping values to the 16-bit range if required.
 *
 * @details Integer variant of hp_fits_write_subset() for the mask output image.
 * Calls hp_fits_correct_data_int() then writes via fits_write_subset_int().
 * Parameters are identical to the float version except @p data is an int array.
 *
 * @param fptr      Open cfitsio fitsfile pointer (write mode).
 * @param group     Group parameter (0 for primary image).
 * @param naxis     Number of image axes.
 * @param naxes     Array of axis sizes.
 * @param data      Full rPixX×rPixY int image; values clamped in place.
 * @param status    cfitsio status code; passed through.
 * @param makeShort If non-zero, clamp to 16-bit short range.
 * @param bZero     FITS BZERO.
 * @param bScale    FITS BSCALE.
 * @param fpixelX   First pixel column to write (1-indexed FITS convention).
 * @param fpixelY   First pixel row to write.
 * @param lpixelX   Last pixel column to write.
 * @param lpixelY   Last pixel row to write.
 * @param xArrayLo  Column offset within the @p data row for the first pixel.
 * @param yArrayLo  Row offset within @p data for the first row.
 * @return cfitsio status code.
 */
int hp_fits_write_subset_int(fitsfile *fptr, long group, long naxis, long *naxes,
                             int *data, int *status, int makeShort,
                             float bZero, float bScale,
                             int fpixelX, int fpixelY, int lpixelX, int lpixelY, int xArrayLo, int yArrayLo) {
    
    int   pixX, pixY, y, x;
    int  *dptr;
    long  fpixel[2], lpixel[2];
    
    hp_fits_correct_data_int(data, (rPixX*rPixY), bZero, bScale, makeShort);
    
    pixX = lpixelX - fpixelX + 1;
    pixY = lpixelY - fpixelY + 1;
    
    fpixel[0] = fpixelX;
    lpixel[0] = lpixelX;
    
    dptr = data;
    /* get to the first row */
    for (y = 0; y < yArrayLo; y++)
        for (x = 0; x < rPixX; x++, dptr++);
    
    for (y = 0; y < pixY; y++) {
        
        fpixel[1] = fpixelY + y;
        lpixel[1] = fpixel[1];
        
        /* get to first pixel in row to write */
        for (x = 0; x < xArrayLo; x++, dptr++);
        
        /* this works! */
        fits_write_subset_int(fptr, group, naxis, naxes, fpixel, lpixel, dptr, status);
        
        /* clear the row, as fits_write_subset does not increment dptr */
        for (x = xArrayLo; x < rPixX; x++, dptr++);
    }
    return *status;
}

/**
 * @brief Fill a float pixel array with a constant value.
 *
 * @param data   Float array to be filled; length must be at least nPixX*nPixY.
 * @param value  Fill value.
 * @param nPixX  Width (used together with nPixY to determine array length).
 * @param nPixY  Height.
 */
void fset(float *data, double value, int nPixX, int nPixY) {
    int    pixelIdx;
    float *pixelPtr;
    pixelPtr = data;
    for (pixelIdx = nPixX*nPixY; pixelIdx--; )
        *(pixelPtr++) = value;
}

/**
 * @brief Fill a double pixel array with a constant value.
 *
 * @param data   Double array to be filled; length must be at least nPixX*nPixY.
 * @param value  Fill value.
 * @param nPixX  Width (used together with nPixY to determine array length).
 * @param nPixY  Height.
 */
void dfset(double *data, double value, int nPixX, int nPixY) {
    int     pixelIdx;
    double *pixelPtr;
    pixelPtr = data;
    for (pixelIdx = nPixX*nPixY; pixelIdx--; )
        *(pixelPtr++) = value;
}

/**
 * @brief Print a cfitsio error report to stderr and terminate the program.
 *
 * @details Calls fits_report_error() to print the full cfitsio error stack and
 * then exits with @p status as the return code.  A no-op if status == 0.
 *
 * @param status  cfitsio status code; non-zero values trigger the error report
 *                and program exit.
 */
void printError(int status) {
    if (status) {
        fits_report_error(stderr, status); /* print error report */
        exit( status );    /* terminate the program, returning error status */
    }
    return;
}

#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/**
 * @brief Portable uniform pseudo-random number generator returning values in
 *        (0, 1).
 *
 * @details Implements the combined multiplicative linear congruential generator
 * described in Numerical Recipes (Press et al.), using three independent
 * congruential sequences shuffled via a Bays–Durham table of 97 entries to
 * remove low-order serial correlations.  The generator is seeded by passing
 * *idum < 0 or on the first call; subsequent calls with *idum > 0 continue the
 * sequence.  Used by getStampStats3() for random pixel sampling when estimating
 * histogram bin widths.
 *
 * @param idum  Seed pointer: pass a negative integer to (re-)seed; *idum is
 *              set to 1 after seeding and must not be altered between calls.
 * @return Pseudo-random double in the open interval (0, 1).
 */
double ran1(int *idum) {
    static long ix1,ix2,ix3;
    static double r[98];
    double temp;
    static int iff=0;
    int j;
    /* void nrerror(char *error_text); */
    
    if (*idum < 0 || iff == 0) {
        iff=1;
        ix1=(IC1-(*idum)) % M1;
        ix1=(IA1*ix1+IC1) % M1;
        ix2=ix1 % M2;
        ix1=(IA1*ix1+IC1) % M1;
        ix3=ix1 % M3;
        for (j=1;j<=97;j++) {
            ix1=(IA1*ix1+IC1) % M1;
            ix2=(IA2*ix2+IC2) % M2;
            r[j]=(ix1+ix2*RM2)*RM1;
        }
        *idum=1;
    }
    ix1=(IA1*ix1+IC1) % M1;
    ix2=(IA2*ix2+IC2) % M2;
    ix3=(IA3*ix3+IC3) % M3;
    j=1 + ((97*ix3)/M3);
    /* if (j > 97 || j < 1) nrerror("RAN1: This cannot happen."); */
    temp=r[j];
    r[j]=(ix1+ix2*RM2)*RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

/**
 * @brief Initialise an index array and sort it by the values in @p list using
 *        an in-place quicksort.
 *
 * @details Fills index[0..n-1] with 0..n-1, then calls quick_sort_1() to
 * sort the index array so that list[index[0]] <= list[index[1]] <= ... <=
 * list[index[n-1]].  The original @p list is not modified.  Used by
 * getPsfCenters() to rank candidate substamp positions by brightness.
 *
 * @param list   Array of n double values to sort (read-only).
 * @param index  Output: index array of length n; initialised and sorted on
 *               return.
 * @param n      Number of elements.
 */
void quick_sort (double *list, int *index, int n) {

    int indexIdx;
    void quick_sort_1();

    for (indexIdx = 0; indexIdx < n; indexIdx++) index [indexIdx] = indexIdx;
    quick_sort_1 (list, index, 0, n-1);

    return;
}



/**
 * @brief Recursive in-place quicksort on an index array, using a median-of-one
 *        pivot.
 *
 * @details Partitions the subarray index[left_end..right_end] around the pivot
 * list[index[mid]] (mid = (left_end+right_end)/2), then recursively sorts each
 * partition.  The sort is ascending in list values.  Called exclusively by
 * quick_sort().
 *
 * @param list       Array of double values (read-only); indexed through @p index.
 * @param index      Index array being sorted in place.
 * @param left_end   Left boundary of the subarray to sort (inclusive).
 * @param right_end  Right boundary of the subarray to sort (inclusive).
 */
void quick_sort_1(double *list, int *index, int left_end, int right_end) {
    int rightIdx,leftIdx,tempIndex;
    double pivotValue;


    pivotValue = list[index[(left_end + right_end)/2]];
    rightIdx = left_end-1;
    leftIdx = right_end+1;

    for(;;) {
        while(list [index[++rightIdx]] < pivotValue);
        while(list [index[--leftIdx]] > pivotValue);
        if (rightIdx < leftIdx){
            tempIndex=index [leftIdx];
            index [leftIdx] = index [rightIdx];
            index [rightIdx] = tempIndex;}
        else if (rightIdx == leftIdx) {
            ++rightIdx; break;}
        else break;
    }

    if (left_end < leftIdx)  quick_sort_1 (list, index, left_end, leftIdx);
    if (rightIdx < right_end) quick_sort_1 (list, index, rightIdx, right_end);

    return;
}

/**
 * @brief Comparator function for qsort() that orders double values in ascending
 *        order.
 *
 * @details Returns -1, 0, or 1 depending on whether *x < *y, *x == *y, or
 * *x > *y respectively.  Passed to qsort() in getStampStats3() to sort a
 * sample of pixel values for interquartile range estimation.
 *
 * @param x  Pointer to the first double value.
 * @param y  Pointer to the second double value.
 * @return -1 if *x < *y, 0 if equal, 1 if *x > *y.
 */
/**** comparison call for the qsort ****/
int flcomp(xVal,yVal)
     double *xVal,*yVal;
{
    if (*xVal>*yVal) return(1);
    else if (*xVal==*yVal) return(0);
    else return(-1);
}

/**
 * @brief Return the smaller of two integers.
 *
 * @param a  First integer.
 * @param b  Second integer.
 * @return The minimum of a and b.
 */
int imin(int a, int b) {
    if (a < b) { return a; }
    else { return b; }
}

/**
 * @brief Return the larger of two integers.
 *
 * @param a  First integer.
 * @param b  Second integer.
 * @return The maximum of a and b.
 */
int imax(int a, int b) {
    if (a < b) { return b; }
    else { return a; }
}

