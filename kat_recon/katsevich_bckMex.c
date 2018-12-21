/*=================================================================
 * katsevich_bckMex.c - Katsevich backprojection in C
 *
 * Function to perform the backprojection of an already filtered
 * projection set using the Katsevich algorithm for cone-beam
 * helicoidal CT.
 *
 * Translation to C of:
 * function katsevich_bckprj(prjFile,outFile,nPrjs,prjSize,pxSize,volSize,volOr,vxSize,angSpan,startAng,P,R0,D)
 *
 * Input:   xxxxxx
 * Output:  xxxxxx
 *
 * 2013 LIM - BiiG - UC3M
 * $Revision: 0.0$
 * $Revision: 0.1$
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "matrix.h"

#ifndef PI
#define PI 3.1415926535f
#endif

#define numArgs 13

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    /* Vars */
    
    /* Loops */
    int l, ix, iy, iz; /* Loop index vars */
    int pos_l, pre_l; /* Loop circular index vars: upper_index, lower_index*/
    /* Pointers to arrays */
    /* Sizes */
    double *prjSize;    /* Projection size */
    double *pxSize;     /* Pixel size */
    double *volSize;    /* Volume size */
    double *vxSize;     /* Voxel size */
    double *volOr;      /* Volume origin */
    /* Angular */
    double *angSpan;    /* Angular positions */
    double *R0;         /* Helix radius */
    float *angCoords;   /* Angular positions in radians */
    /* Storage */
    float *prjStack;
    float *volStack;
    /* Files */
    FILE *fdPrj;       /* Projection data file */
    FILE *fdVol;       /* Volume output file */
    /* Strings */
    char *prjFile;     /* Projection file */
    char *outFile;     /* Output file */
    
    /* Aux vars: miscellaneous */
    int tmp;
    int nPrjs, uSize, wSize, xSize, ySize, zSize, sliSize; /* Number of projections and sizes, i.e. loop limits */
    float  Dsq, pi2;
    float *mag;
    float *TD_cte;
    float g_int; /* Interpolated projection value */
    float x_0, y_0, z_0; /* Offsets */
    float uStep, wStep, xStep, yStep, zStep; /* Steps */
    float angSpan_rd, max_rd, lambdaStart, lambdaEnd; /* Angular stuff */
    float *lambdaStep;
    /* Aux vars for processing inside the loop */
    float lambdaCur, lambdaCos, lambdaSin, lambdaPre, lambdaCosPre, lambdaSinPre, lambdaPos, lambdaCosPos, lambdaSinPos;
    float x, y, z;    /* Volume coords */
    float u_st, v_st, w_st, u_st_pre, v_st_pre, u_st_pos, v_st_pos, w_st_inc; /* Detector indexes */
    float w_bot, w_top, z_bot, z_top, w_bot_pos, w_top_pre, z_bot_pos, z_top_pre; /* Window limits */
    float sca, weigh; /* Scaling factor and weights */
    /* Addressing the detector */
    float u_tmp, w_tmp, u_res, w_res;
    int u_ind, w_ind;
    /* Input data */
    float  startAng, P, D;
    
    
    /* Process */
    /* Check arg number */
    if (nrhs != numArgs)
    {
        mexErrMsgTxt("ERROR: Wrong number of arguments\n");
    }
    /* Extract the args */
    /* Filenames */
    /* Input 1 and 2 must be strings */
    if ((mxIsChar(prhs[0]) != 1) || (mxIsChar(prhs[1]) != 1))
        mexErrMsgTxt("Input must be a string.");
    
    /* Input 1 and 2 must be row vectors */
    if ((mxGetM(prhs[0])!=1) || (mxGetM(prhs[1])!=1))
        mexErrMsgTxt("Input must be a row vector.");
    
    /* Get the strings */
    prjFile = mxArrayToString(prhs[0]);
    outFile = mxArrayToString(prhs[1]);
    
    /* Now the rest of the args */
    /* Scalars */
    nPrjs    = (int) mxGetScalar(prhs[2]);
    /*angSpan  = (float) mxGetScalar(prhs[8]);*/
    startAng = (float) mxGetScalar(prhs[9]);
    P        = (float) mxGetScalar(prhs[10]);
    D        = (float) mxGetScalar(prhs[12]); // AO - SDD considering mean value of the calib file
    /*AO - Transform into array args */
    /* R0       = (float) mxGetScalar(prhs[11]);*/
    
    
    /* Arrays */
    prjSize = mxGetPr(prhs[3]);     // [N_rows, N_cols] (px)
    pxSize  = mxGetPr(prhs[4]);     // [S_row S_col] (mm)
    volSize = mxGetPr(prhs[5]);     // [N_rows, N_cols, N_slices] (vx)
    volOr   = mxGetPr(prhs[6]);     // [or_x or_y or_z]
    vxSize  = mxGetPr(prhs[7]);     // [S_row S_col S_slice] (mm)
    angSpan = mxGetPr(prhs[8]);     // [ang0 ang1 ...angnPrjs-1] (deg)
    R0      = mxGetPr(prhs[11]);    // [SAD0 SAD1 ...SADnPrjs-1] (mm)
    
    /* Projection sizes */
    uSize = (int) prjSize[0];       // N_rows
    wSize = (int) prjSize[1];       // N_cols
    uStep = (float) pxSize[0];      // S_row
    wStep = (float) pxSize[1];      // S_col
    
    /* Volume sizes */
    xSize = (int) volSize[0];
    ySize = (int) volSize[1];
    zSize = (int) volSize[2];
    sliSize = xSize*ySize;
    xStep = (float) vxSize[0];
    yStep = (float) vxSize[1];
    zStep = (float) vxSize[2];
    x_0   = (float) volOr[0]; /* For the time being , we assumed this to be zero, only offset in z */
    y_0   = (float) volOr[1]; /* For the time being , we assumed this to be zero, only offset in z */
    z_0   = (float) volOr[2]; /* Offset in z */
    
    mexPrintf("Prj Size: %dx%d\n",uSize,wSize);
    mexPrintf("Vol Size: %dx%dx%d\n",xSize,ySize,zSize);
    mexPrintf("Vol orig: %.2fx%.2fx%.2f\n",x_0,y_0,z_0);
    mexEvalString("drawnow");
    
    /* AO - Substitute the angular stimation by the vector of angles from the acquisition IN RAD
     * angCoords[] = angles*PI/180.0f
     */
    /* Angular params - Be careful with deg-rad mistakes */
//     angSpan_rd  = (angSpan.*PI)./180.0f;
//     lambdaStep  = angSpan_rd/nPrjs; /* AO - mean from diff between prev and post angle? therefor lambdaStep[]*/
//     lambdaStart = (startAng*PI)/180.0f;
//     lambdaEnd   = lambdaStart + angSpan_rd - lambdaStep;
//     for (l = 0; l < nPrjs; l++)
//     {
//         angCoords[l] = lambdaStart + l*lambdaStep;
//     }
    /* Generate vector for angular coords */
    angCoords = (float*) mxCalloc(nPrjs,sizeof(float));
    lambdaStep = (float*) mxCalloc(nPrjs,sizeof(float));
    for (l = 0; l < nPrjs; l++)
    {
        angCoords[l] = (float)angSpan[l];
        angCoords[l] = (float)((angCoords[l]*PI)/180.0f);
    }
    /* Constant quantities */
    TD_cte = (float*) mxCalloc(nPrjs,sizeof(float));
    mag = (float*) mxCalloc(nPrjs,sizeof(float));
    
    max_rd = (2*PI);
    for (l = 0; l < nPrjs; l++)
    {
        pos_l = l + 1 >= nPrjs ? (l+1) - nPrjs : l+1;
        pre_l = l-1 < 0 ? (l-1) + nPrjs : l-1;
        if(l == 0)
        {
            lambdaStep[l]=fabsf((angCoords[pos_l] - angCoords[l])+(angCoords[l] + (max_rd-angCoords[pre_l])))/2.0f;
        }else if(l == nPrjs-1)
        {
            lambdaStep[l]=fabsf(((max_rd+angCoords[pos_l]) - angCoords[l])+(angCoords[l] - angCoords[pre_l]))/2.0f;
            
        }else
        {
            lambdaStep[l]=fabsf((angCoords[pos_l] - angCoords[l])+(angCoords[l] - angCoords[pre_l]))/2.0f;
        }
//         if(l == nPrjs-1)
//             lambdaStep[l]=fabsf(((max_rd+angCoords[pos_l]) - angCoords[l]));
//         else
//             lambdaStep[l]=fabsf((angCoords[pos_l] - angCoords[l]));
        TD_cte[l] = (float)(P/((2*PI*D)*((float)R0[l])));
        // AO - mag is an array
        mag[l]    = (float) (D/(float)R0[l]); // AO - variable magnification due to R0 variability
//          mexPrintf("ang span %f ang coord %f lambda step %f\n",angSpan[l],angCoords[l], lambdaStep[l]);
//         mexPrintf("TD_cte %f mag %f\n",TD_cte[l], mag[l]);
//          mexEvalString("drawnow");
//         mexPrintf("Data %d of %d \t step %f\n",l,nPrjs, lambdaStep[l]);
//         mexEvalString("drawnow");
    }
    
    Dsq    = D*D; /*AO- distance square law*/
    pi2    = PI/2.0f;
    
    
    /* Storage */
    prjStack = (float*) mxCalloc(uSize*wSize,sizeof(float)); /* We need at least three projections to allow the backprojection according to Noo et al. */
    volStack = (float*) mxCalloc(sliSize*zSize,sizeof(float));
    
    /* Open files */
    if((fdPrj=fopen(prjFile,"rb"))==NULL)
    {
        mexErrMsgTxt("ERROR: cannot open projections file\n");
    }
    
    /* Loop through projections, then, for each projection through the entire volume */
    
    for (l = 0; l < nPrjs; l++)
//      for (l = 0; l < 1; l++)
    {
        
        /* Report progress */
        mexPrintf("Processing projection %d of %d\n",l,nPrjs);
//      Debugging
//         mexPrintf("ang coord %f lambda step %f\n",angCoords[l], lambdaStep[l]);
//         mexPrintf("TD_cte %f mag %f\n",TD_cte[l], mag[l]);
        mexEvalString("drawnow");
//      End of debug
        /* Read data and update index */
        tmp = (int) fread(prjStack,uSize*wSize,sizeof(float),fdPrj);
        if(tmp<=0)
        {
            mexErrMsgTxt("ERROR: read %d\n", tmp);
        }
        /* Calculations dependent only on angular position */
        lambdaCur = angCoords[l];
        lambdaCos = cosf(lambdaCur);
        lambdaSin = sinf(lambdaCur);
        
        /* Previous and next angles */
        pos_l = l + 1 >= nPrjs ? (l+1) - nPrjs : l+1;
        pre_l = l-1 < 0 ? (l-1) + nPrjs : l-1;
//         mexPrintf("Post pos %d pre pos %d\n", pos_l, pre_l);
//         mexEvalString("drawnow");
        /* AO - modify to get previous and next from angCoords[l-+1] CIRCULAR VECTOR: previous angle for 0º is last angles value*/
//         lambdaPre    = angCoords[l] - lambdaStep[l];/* AO - angCoords[l-1]*/
        lambdaPre    = angCoords[pre_l];
        lambdaCosPre = cosf(lambdaPre);
        lambdaSinPre = sinf(lambdaPre);
//         lambdaPos    = angCoords[l] + lambdaStep[l]; /* AO - angCoords[l+1]*/
        lambdaPos    = angCoords[pos_l];
        lambdaCosPos = cosf(lambdaPos);
        lambdaSinPos = sinf(lambdaPos);
        
//      Debugging
//         mexPrintf("LambdaCur %d LamdaPos %d lambdaPre %d \n ",lambdaCur,lambdaPos,lambdaPre);
//         mexEvalString("drawnow");
//      end of debugging
        
        /* Volume loop - first x, then y, last z -> See Noo etal for the reason */
        for (ix = 0; ix < xSize; ix++)
        {
            /* Calculations not influenced by y neither z */
            x = (ix - xSize/2 - 1) * xStep; /* x coord */
            
            /* Y loop */
            for (iy = 0; iy < ySize; iy++)
            {
                /* Calculations not influenced by z */
                y = (iy - ySize/2 - 1) * yStep; /* y coord */
                /* Detector indexes */
                /* AO - R0[l]*/
//                 v_st = R0 - x*lambdaCos - y*lambdaSin; /* st stands for star, keep notation coherent with refs */
                v_st = (float)(R0[l] - x*lambdaCos - y*lambdaSin); /* st stands for star, keep notation coherent with refs */
                u_st = (D/v_st)*(-x*lambdaSin + y*lambdaCos);
                /* Also for previous and post-projection */
//                 v_st_pre = R0 - x*lambdaCosPre - y*lambdaSinPre;  /* AO - R0[l-1]*/
                v_st_pre = (float)(R0[pre_l] - x*lambdaCosPre - y*lambdaSinPre);
                u_st_pre = (D/v_st_pre)*(-x*lambdaSinPre + y*lambdaCosPre);
//                 v_st_pos = R0 - x*lambdaCosPos - y*lambdaSinPos; /* AO - R0[l+1]*/
                v_st_pos = (float)(R0[pos_l] - x*lambdaCosPos - y*lambdaSinPos);
                u_st_pos = (D/v_st_pos)*(-x*lambdaSinPos + y*lambdaCosPos);
                
                /* Tam-Danielsson window limits */
                w_bot = -TD_cte[l] * (u_st*u_st + Dsq) * (pi2 + atanf(u_st/D));
                w_top =  TD_cte[l] * (u_st*u_st + Dsq) * (pi2 - atanf(u_st/D));
                z_bot = (v_st/D)*w_bot + (P/(2*PI))*lambdaCur;
                z_top = (v_st/D)*w_top + (P/(2*PI))*lambdaCur;
                
                /* Compute limits for previous and post-projections and convert them to z */
                w_bot_pos = -TD_cte[pos_l] * (u_st_pos*u_st_pos + Dsq) * (pi2 + atanf(u_st_pos/D));
                w_top_pre =  TD_cte[pre_l] * (u_st_pre*u_st_pre + Dsq) * (pi2 - atanf(u_st_pre/D));
                
                /* Convert to z */
                z_bot_pos = (v_st_pos/D)*w_bot_pos + (P/(2*PI))*lambdaPos; /*AO - angCoords[l+1]*/
                z_top_pre = (v_st_pre/D)*w_top_pre + (P/(2*PI))*lambdaPre; /*AO - angCoords[l-1]*/
                
                /* Compute constant scale factor */
//                 sca = lambdaStep/(2*PI*v_st); /* AO - lambdaStep[l]*/
                sca = (float)(lambdaStep[l]/(2*PI*v_st));
                /* Z loop - For this we use the cone beam cover method described */
                /* by Fontaine */
                /* Not sure it is so accurate as the PI-line segment approach */
                /* but it is much more efficient */
//                 z        = (zSize/(mag))*zStep - z_0; /* AO - mag[l]*/
                z        = (float)(zSize/(mag[l]))*zStep - z_0;
                w_st     = (D/v_st)*(z - (P/(2*PI))*lambdaCur);
                w_st_inc = (D/v_st)*zStep;
                
                /* Speedup the loop, loop only through the valid z */
                for (iz = 0; iz < zSize; iz++)
                {
                    z = z - zStep; /* z coord */
                    /* Detector index */
                    w_st = w_st - w_st_inc;
                    
                    /* Is the pixel within the Tam-Danielsson window? */
                    /* or within the PI-line "extended" limits */
                    if ((w_st > w_bot) && (w_st < w_top))
                    {
                        /* Get weight */
                        if (z < z_bot_pos)
                            weigh = (float) (fabs(z_bot - z)/fabs(z_bot - z_bot_pos) + 0.5f);
                        else if (z > z_top_pre)
                            weigh = (float) (fabs(z_top - z)/fabs(z_top - z_top_pre) + 0.5f);
                        else
                            weigh = 1.0f;
                        
                        /* Get the interpolated value */
                        w_tmp = w_st/wStep + wSize/2;
                        u_tmp = u_st/uStep + uSize/2;
                        w_ind = (int)(w_tmp);
                        u_ind = (int)(u_tmp);
                        w_res = w_tmp - w_ind;
                        u_res = u_tmp - u_ind;
                        
                        
                        /* If outside of projection, get out of here */
                        if ((u_ind < 0) || (u_ind >= uSize - 1) || (w_ind < 0) || (w_ind >= wSize - 1))
                            continue;
                        
                        g_int = prjStack[w_ind*uSize + u_ind]*(1.0f-u_res)*(1.0f-w_res) + prjStack[w_ind*uSize + (u_ind+1)]*u_res*(1.0f-w_res) + prjStack[(w_ind+1)*uSize + u_ind]*(1.0f-u_res)*w_res + prjStack[(w_ind+1)*uSize + (u_ind+1)]*u_res*w_res;
                        
                        /* Update voxel */
                        volStack[iz*sliSize + iy*xSize + ix] = volStack[iz*sliSize + iy*xSize + ix] + sca*weigh*g_int;
                        
                    } /* End if */
                } /* End for z */
            } /* End for y */
        } /* End for x */
        /* Move to next prj */
    } /* End for l */
    
    /* Close prj file */
    fclose(fdPrj);
    
    /* Save volume file */
    fdVol = fopen(outFile,"wb");
    fwrite(volStack,sliSize*zSize,sizeof(float),fdVol);
    fclose(fdVol);
    
    /* Free resources */
    mxFree(angCoords);
    mxFree(prjStack);
    mxFree(volStack);
    mxFree(prjFile);
    mxFree(outFile);
    
}