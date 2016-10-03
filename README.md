# gencode
#commands to make this stuff
gcc -g -o gencode gencode.c submalloc.c
gcc -g -o gencode -DNOSUBMALLOC gencode.c submalloc.c
gcc -g -o gencode -DSNIPLISTPT gencode.c submalloc.c
gcc -g -o tstsubmal tstsubmal.c submalloc.c
gcc -g -o tstsubmal.nosub -DNOSUBMALLOC tstsubmal.c submalloc.c

gencode is a program that can reconstruct long genetic codes from
snippets.  The snippets can't have errors, and the original code
can't have long repeats.

gencode has two modes, one mode to split multiple copies of a long
generated random sequence into snippets and the other to put the
snippets back together.  The snippets are generated in random order
and in one of 4 versions.  The versions are forward, forward mirror
helix, reverse, and reverse mirror helix.  The solution is output in
those 4 versions too.

4 copies of the code in 20-100 letter snippets can be solved with
minolap of 10 for codes up to about 100K long.  5 copies of the code
in 100-400 letter snippets can be solved with minolap of 10 for codes
up to about 10M long.  The 10M case works better with a minolap of 12. 
6 copies of the code in 100-400 letter snippets with minolap of 16 is
needed for codes 100M long.  The program can handle a case of 6 
copies of a 650M letter code in 100-400 long snippets on a 16G linux
machine.

Type the command without arguments to get the command syntax.  To
generate sinppets, type:

	./gencode 0 filebase -l codelen -i minsniplen -x maxsniplen -c ncopy

It will generate filebase.code, which has the 4 versions of the
original code, and filebase.snip, which has all the snippets. 
Filebase.snip has a line for each snippet with the length and snippet
string.  Each snippet string will be in one of the 4 versions.  If you
add a "-?" argument to the command, it will find the longest duplicate
string of letters in the code.

To reconstruct the original code form the snippets, type:

	./gencode 1 filebase -o minolap -d 4 >filebase.sol

It will look in filebase.snip for snippets.  The default minolap is
10.  It should be increased to 12 for codes longer than 10 million
letters, 14 for codes longer than 100 million, and 16 for codes longer
than 300 million.  The "-d 4" argument prints some memory statistics. 
Gencode first reads the snippets and makes the 4 version copies of
each snippet.  It then eliminates snippets that are wholely contained
in other snippets.  It then works by making several passes of finding
the best new overlap fit for each snippet.  It combines snippets into
ropes of snippets.  It adds to ropes and combines ropes.  At then end
it may have a bunch of ropes that could not be joined because there is
insufficient overlap between one rope and the next.  This is because
all the copies were cut too close to the same point in the code.  If
there are not too many ropes, it goes through an end join process that
tries to pick the best rope sequences to join.
