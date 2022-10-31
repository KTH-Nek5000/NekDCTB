## Lossless compression with BZIP2

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

-   **Run-length encoding.** Where portions of data are represented as a
    single data value and a count instead of the original, for example
    representing a string of \"5555222\" as \"4532\". This is
    particularly useful for the truncated velocity fields as many
    entries has been set to zero.

-   **Burrows--Wheeler transform**. Which is a reversible transformation
    that rearranges the entries of a character string into runs of
    repeated characters. For example, a string \"\*BANANA\|\" will be
    rearranged as \"BNN\*AA\|A\", where more identical characters can be
    seen next to each other.

-   **Move-to-front transform**. The main idea of the algorithm is that
    each symbol is replaced by its index in a list of most used symbols,
    if a symbol appears in sequence it would be replaced by zeros while
    if it does not appear frequently is replaced by a large number.
    Applying the algorithm to the string \"bananaaa\" having as a list
    of symbols the letters in the alphabet, will produce the following
    sequence: \"1,1,13,1,1,1,0,0\" .

-   **Run-length encoding.** once again.

-   **Huffman coding.** The algorithm produces a code table for encoding
    the source symbol. The table is obtained by calculating the
    probability (or frequency of occurrence) for each value in the
    source to occur. Symbols with high probability will be represented
    with codes that use less bits that those used for less common
    symbols.
