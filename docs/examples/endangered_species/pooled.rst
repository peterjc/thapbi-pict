Universal animal and plant DNA barcodes and mini-barcodes
=========================================================

We have very briefly reviewed the output of each of the animal and plant
markers, noting many have no sequences at the THAPBI PICT default minimum
abundance threshold. Now we discuss the pooled results produced by the
``run.sh`` shell script (which literally pooled the markers for each sample
by concatenating the intermediate files together).

Please open the ``pooled.samples.onebp.xlsx`` sample report, and to look at
the unknown reads see ``pooled.reads.onebp.xlsx``. Overall the replicates
are reassuringly consistent.

The automated model assessment output in ``pooled.assess.onebp.tsv`` is
also worth review. Note this only looks at the experimental mixtures where
there is a ground truth - not the traditional medicine samples where the
true species content is unknown (detection of human and rat DNA included).

There are false positives for *Acipenser schrenckii*, human, and alternative
genus level matches in *Brassica* and *Lactuca* (as discussed in the paper).

If the sample database had been more inclusive there would have been many
more false positives. For example, the trnL-UAA sequence perfectly matching
AP007232.1 *Lactuca sativa* is also a perfect match for MK064549.1 *Luisia
teres*. Similarly, the Mini-rbcL sequence perfectly matching AP012989.1
*Brassica nigra* and MG872827.1 *Brassica juncea* also matches MN056359.2
*Raphanus sativus* (and more). This demonstrates the difficulties in curating
an appropriate marker database - and the content should depend in part on your
target samples.

Currently the provided references sequences (and thus classification databases
used) lack any markers for *Beta vulgaris*, *Cycas revoluta*, *Parapenaeopsis*
sp., *Triticum aestivum* or *Zea mays*. Most of these were present at only a
few percent dry weight, and are likely present below the default minimum
abundance threshold. This explains the false negatives.

It appears that the THAPBI PICT default minimum abundance threshold of 100
reads is too stringent for detecting all the markers in a complex pool like
this. Including negative sequencing controls would help set an objective
lower bound.
