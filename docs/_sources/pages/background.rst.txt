============================================
Entropy and information theory
============================================

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

==================================================================
Steps for data compression in turbulent fields.
==================================================================

Given that turbulence is a multi-scale and seemlingly random phenomenon,
data sets that describe it tend to have quite a large entropy. Making
lossless compression ineficient. It is for this reason that in order to
achieve high compression ratios, one alternatie is to first include
steps that reduce the entropy, generally converting the compression into
a lossy process.

Lossy data truncation
------------------------------------


There are many ways to perform lossy compression, however,
transfromation-based aproaches are among the most used. JPEG for image
compression, for example, uses this aproach. Given that the data in
Nek5000 is distributed in a mesh with certain characteristics, the
aproach taken in this library was to truncate the data in the legendre
space while controlling the error. This is benefitial given that the
user has total control in the fidelity of the data after it is stored.

Legendre Polynomials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spectral element solvers such as Nek5000, Neko and NekRS by discretize the domains with high
order finite elements, each of them containing a number of
Gauss-Lobatto-Legendre collocation points that allow them to perform
numerical differentiation and integration thanks to numerical
quadrature. A good characteristic of these points is that a set of
polynomials known as the **Legendre polynomials** form
an orthogonal basis :math:`\phi` around them, therefore the value
:math:`u` of a known function at point :math:`x` can be represented as a
sum of the values of the orthogonal basis at different orders scaled by
a set of Legendre coefficients :math:`\hat{u}` as follow:

.. math::

       
   u(x)=\sum_{i=0}^{N}\hat{u}_i\phi_i(x)

The polynomials of order 0 and 1 are :math:`\phi_0(x)=1` and
:math:`\phi_1(x)=x`, respectively, and the rest are obtained with the
recursive formula:

.. math::


       (i+1)\phi_{i+1}(x)=(2i+1)\phi_{i}(x)-(i)\phi_{i-1}(x), \quad i=1,..N,

For convenience we also define now the integration weight matrix
produced by the GLL quadrature:

.. math::


       W_{ij}=\frac{2}{N(N+1)}\frac{1}{\phi^2_N(x_i)}\delta{ij}.

Although this equation works for the value in a single
point, in numerical methods we are often interested in a set of them.
Therefore, if we define :math:`\textbf{u}=[u_0, u_1, ... u_{N-1}]`, it
is possible to transform the entire array by using a transformation
matrix as follows:

.. math::


       \textbf{u}=T\hat{\textbf{u}},

where :math:`T_{ij}^T=f_j\phi_{i}(x_j)` and :math:`f_j` are scaling
factors that make :math:`T` orthogonal with respect to the integration
weights such that :math:`T^TWT=I` and are defined as:

.. math::


   f_j=
    \begin{cases} 
         \sqrt{(2(j+1)-1)/2} & 0\leq j\leq N-1 \\
         \sqrt{(N-1)/2} & j=N 
    \end{cases}

Based on the previous explanation, it is straight forward to realize
that a physical signal :math:`u` has an equivalent representation in the
Legendre space :math:`\hat{u}=T^{-1}u`, the main difference however, is
that not all entries in :math:`\hat{u}` have the same relevance in the
physical representation of the signal. Therefore, discarding or
approximating some of the Legendre coefficients will provide fairly
accurate physical reconstruction of the signal as long as some key
elements are maintained. It is precisely this characteristic that makes
this sort of approach a relevant case study for compression algorithms
which is what we will focus in the following sections. We will initially
focus in the effect that discarding information has in the reconstructed
signal, and then analyze what needs to be stored in order for the
algorithms to effectively compress.

Down-sampling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In signal processing, often times the high frequencies can not be
perceived by humans or are not relevant for the signal representation in
physical space. Turbulence is not an exception to this characteristic,
as it is a process with high dissipation primarily at the smallest
eddies, i.e, the higher frequencies. For this reason a natural choice of
information to discard is is precisely that contained in higher
frequencies.

As mentioned before, in Legendre space, the information about the signal
is contained mainly in the coefficients :math:`\hat{u}`. The process of
down-sampling consist of filtering a number of the coefficients that
represent the importance of the high order Legendre polynomials (the
higher frequencies) and transforming back such that an operation as the
following is performed:

.. math::


       \tilde{u}=T^{-1}FTu,

Where :math:`F` is a filtering operator that consist of a "broken
identity matrix". Eq. `[eq:ds-1] <#eq:ds-1>`__ summarize 3 main steps:

1. Transform the signal into the Legendre space: :math:`\hat{u}=Tu`.

2. Filter the signal: :math:`\tilde{\hat{u}}=F\hat{u}`.

3. Transform the down-sampled coefficients back to physical space:
   :math:`\tilde{u}=T^{-1}\tilde{\hat{u}}`.

It is important to note 2 main points: first that
:math:`\tilde{u} \neq u` unless :math:`F=I`, so down-sampling
undoubtedly brings with it some errors and second that the down-sampled
coefficients contain less information than the original ones since only
:math:`k` coefficients are kept, while the rest are set to 0 trhough the
filtering matrix :math:`F`. The content of :math:`\tilde{\hat{u}}` will
be an array with entries arrayed in a similar fashion as shown now:

.. math::


   \tilde{\hat{u}}=
   \begin{bmatrix}
   \hat{u}_0 \\
   \hat{u}_1 \\
   ...\\
   \hat{u}_{k-1} \\
   \hat{u}_k \\
   0 \\
   .. \\
   0 \\
   \end{bmatrix}.

It is possible to see that since :math:`\tilde{\hat{u}}` has many zeros
in it, it is much easier to represent it in memory than the
reconstructed signal :math:`\tilde{u}` and through lossless compression
algorithms that include run-length encoding or Huffman encoding, the
array in disk will only occupy the size that is needed for the :math:`k`
entries kept while the extra zeros will represent only a minor storage
overhead.

Discreet Legendre truncation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

[@Otero2018] recognized the use of down-sampling but extended it for
turbulence in spectral element methods. They realized that although
higher frequencies should have the lower energies (lower magnitude
coefficients), numerically in the Legendre space, this is not the case,
particularly when considering 3D simulations with the data structures
used in spectral element solvers.

Therefore they introduced extra steps to the algorithm that ensure that
it is not the higher frequencies that are filtered out, instead it is
the frequencies that posses the lower energies. The principles are still
the same as in down-sampling but instead immediately filtering out the
coefficients, first a sorting operation is performed. The algorithm is
as follow:

1. Transform the signal into the Legendre space: :math:`\hat{u}=Tu`.

2. Sort the coefficients :math:`\hat{u}^{s}=sort(\hat{u})`

3. Filter the sorted signal: :math:`\tilde{\hat{u}}^{s}=F\hat{u}^{s}`.

4. Unsort the signal:
   :math:`\tilde{\hat{u}}=unsort(\tilde{\hat{u}}^{s})`.

5. Transform the truncated coefficients back to physical space:
   :math:`\tilde{u}=T^{-1}\tilde{\hat{u}}`.

For this case, the array :math:`\tilde{\hat{u}}` will not possess all
the zeros grouped at the end as would be the case in down-sampling,
however with the correct encoding technique, the storage overhead is
negligible. This method has the advantage that if the same number
:math:`k` of frequencies are kept, lower errors will be obtained due to
the fact that in this method the most influential frequencies are
stored.

Additionally to this procedure, [@Otero2018] also introduced a way to
calculate the reconstruction error at run-time without the need to
transforms back to physical space by taking advantage of 1. how the L2
norm of a physical variable is calculated and 2. the Legendre quadrature
for integration, as follows:

.. math::


       \lVert  u\rVert_{2}=\sqrt{\int_{\Omega}u^2 d\Omega}=\sqrt{u^TWu}=\sqrt{(T\hat{u})^TW(T\hat{u})}=\sqrt{\hat{u}^TW\hat{u}}.

This provided an efficient way to control error an compress only up to a
user predefined threshold :math:`\epsilon` of the error.

The DLT method was initially proposed in a theoretical manner but has
been since extended by to include run-time encoding and produce a fully
functioning lossy compression algorithm that in fact matches well with
the theoretical expectations one would have of such a lossy compression.

We have further extended the code to also control the Peak Signal To
Noise ratio:

.. math::


   PSNR=20 \log{\frac{max(u)-min(u)}{\lVert  u\rVert_{2}}}

Controling this error as well, allows the procedure to control in a more
localized manner where more compression is allowed, in such a way that
in regions with high fluctuations, where the spectrum is flatter, less
error is allowed to mantain fidelity over compression.

Lossless data encoding
------------------------------------

After the entropy of the data set has been reduced by a lossy truncation
step. It can be losslessly compressed without introducing new errors by
a multitude of schemes that serve this purpose. In this project we have
chosen to use bzip2 trough ADIOS2.

For the lossless compression step, we use the library ADIOS2
[@godoy2020adios] in order to have a data management framework that
would allow us to transfer data in an efficient way, apart from
seamlessly integrating many lossless compression algorithms, including
**bzip2**, which we decided to use. In principle, any other form of
lossless compression would serve the same purpose.

The input for this step is a velocity or pressure array that has been
truncated, i.e., many of its entries has been set to be zero. Then the
bzip algorithm applies several layers of compression techniques in the
following order:

-  **Run-length encoding.** Where portions of data are represented as a
   single data value and a count instead of the original, for example
   representing a string of "5555222" as "4532". This is particularly
   useful for the truncated velocity fields as many entries has been set
   to zero.

-  **Burrows–Wheeler transform**. Which is a reversible transformation
   that rearranges the entries of a character string into runs of
   repeated characters. For example, a string ``"*BANANA|"`` will be
   rearranged as "BNN*AA|A", where more identical characters can be seen
   next to each other.

-  **Move-to-front transform**. The main idea of the algorithm is that
   each symbol is replaced by its index in a list of most used symbols,
   if a symbol appears in sequence it would be replaced by zeros while
   if it does not appear frequently is replaced by a large number.
   Applying the algorithm to the string "bananaaa" having as a list of
   symbols the letters in the alphabet, will produce the following
   sequence: "1,1,13,1,1,1,0,0" .

-  **Run-length encoding.** once again.

-  **Huffman coding.** The algorithm produces a code table for encoding
   the source symbol. The table is obtained by calculating the
   probability (or frequency of occurrence) for each value in the source
   to occur. Symbols with high probability will be represented with
   codes that use less bits that those used for less common symbols.


