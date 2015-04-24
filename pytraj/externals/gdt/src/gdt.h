unsigned short * gdtCPU(double *arr, int conformers,int protlen, int score);
unsigned short * gdtCPUOneReference(double *reference, double *arr,  int conformers,int protlen,int score);
unsigned short * gdtCPUOneReferenceExt(double *reference, double *arr,  int conformers,int protlen,int score);
    /*
     *  score:
     * 1 = GDT
     * 2 = TM score
     * 3 = MaxSub
    */
