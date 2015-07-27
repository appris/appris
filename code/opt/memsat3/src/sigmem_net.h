#define IPERGRP (21)

#define WINL (-9)
#define WINR (9)

#define NUM_IN	((WINR-WINL+1)*IPERGRP)	/* number of input units */
#define NUM_HID (15)			/* number of hidden units */
#define NUM_OUT (4) 			/* number of output units */

#define TOTAL		(NUM_IN + NUM_HID + NUM_OUT)

#define IALPHA		(0.9)		/* Initial smoothing factor */
#define ILRATE		(0.001)         /* Initial learning rate */
#define RWMAX		(0.03)		/* Max random weight value */
#define SAMPINTV	(1)		/* Sampling interval */
#define SAMPSIZE	(0)		/* Samples per weight update */

#define CPU_MAX		(10*3600.0)	/* Maximum CPU time (secs) for learning */

#define noEPOCHMODE

