unsigned short * gdtCPU(float *arr, int conformers,int protlen, int score);
unsigned short * gdtCPUOneReference(float *reference, float *arr,  int conformers,int protlen,int score);
unsigned short * gdtCPUOneReferenceExt(float *reference, float *arr,  int conformers,int protlen,int score);
    /*
     *  score:
     * 1 = GDT
     * 2 = TM score
     * 3 = MaxSub
    */
