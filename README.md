A simple Shiny app to deconvolute luciferase reporter data.

Implements the idea from [Brown 2008 - Inferring gene expression dynamics from reporter protein levels](https://pubmed.ncbi.nlm.nih.gov/18956349/) - Biotechnol J
. 2008 Nov;3(11):1437-48. [doi: 10.1002/biot.200800133](https://doi.org/10.1002/biot.200800133).

Can also calculates initial rate, time to maximum and decay rate of the deconvoluted trace.

# Usage

Requires a CSV file as input formatted with two columns called "Time" and "Luminescence"

Input the halflife of your reporter in the same time units as time.

Press "Deconvolute". You're done!

