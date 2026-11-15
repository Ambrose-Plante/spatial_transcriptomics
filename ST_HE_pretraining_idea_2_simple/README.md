# ST_HE_pretraining_idea_2_simple

## Overview

This notebook generalizes the previous H&E-only pretraining idea into a paired histology/spatial-transcriptomics representation model. One encoder compresses H&E imagery while another compresses the very high-dimensional gene-expression tile into a smaller latent feature space. The latent representations can be concatenated or compared before being passed to H&E and ST decoders. Reconstruction losses train each modality to preserve its own information, while cosine-similarity-style ideas in the notebook explore whether the two latent spaces can also be aligned. The intended result is a joint representation in which tissue morphology and local gene expression provide complementary information.

## Algorithmic workflow

1. Encode H&E
2. Encode gene-expression tile
3. Combine/align latent features
4. Reconstruct modalities
5. Calculate reconstruction/alignment losses
6. Optimize jointly
