Theoretical Background
======================

Entropy in information theory
-----------------------------

Based on `Li et
al. <https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13336>`__ in
information theory, the shannon entropy :math:`H` is an optimal lower
bound on the expected number of bits needed to represent symbols
:math:`a_i` in a data set :math:`S`. And it is calculated as:

.. math::


   H=-\sum_{i_1}^n p(a_i)\log_{2}{p(a_i)}

Where :math:`p(a_i)` is the probability of symbol :math:`a_i` to appear
in the data set :math:`S` and the quantity :math:`\log_{2}{p(a_i)}` can
be interpreted as the amount of information that symbol :math:`a_i`
posess. Analyzing the expression, it is possible to observe that
:math:`H` is maximized if all symbols in the data set :math:`S` have the
same probability to appear, and minimized if the probaility of a single
event is 1, while the rest is zero. It is therefore possible to observe
that the entropy is also a measure of the variance or uniformity of the
data set :math:`S`.

A family of compression algorithms exist that try to represent data sets
in such a way that the information is stored at this minimun size
bounded by :math:`H`, however it is clear that if the variance in the
data set is high, such that :math:`H` is high itself, then the size of
data set :math:`S` after being optimally written to the disk might still
be large.

Steps for data compression in turbulent fields.
-----------------------------------------------

Given that turbulence is a multi-scale and seemlingly random phenomenon,
data sets that describe it tend to have quite a large entropy. Making
lossless compression ineficient. It is for this reason that in order to
achieve high compression ratios, one alternatie is to first include
steps that reduce the entropy, generally converting the compression into
a lossy process.

Lossy data truncation
~~~~~~~~~~~~~~~~~~~~~

There are many ways to perform lossy compression, however,
transfromation-based aproaches are among the most used. JPEG for image
compression, for example, uses this aproach. Given that the data in
Nek5000 is distributed in a mesh with certain characteristics, the
aproach taken in this library was to truncate the data in the legendre
space while controlling the error. This is benefitial given that the
user has total control in the fidelity of the data after it is stored.

Details on the theory for this method can be found `following this
link <./truncation.rst>`__

Lossless data encoding
~~~~~~~~~~~~~~~~~~~~~~

After the entropy of the data set has been reduced by a lossy truncation
step. It can be losslessly compressed without introducing new errors by
a multitude of schemes that serve this purpose. In this project we have
chosen to use bzip2 trough ADIOS2. More details can be found `following
this link <./encoding.rst>`__
