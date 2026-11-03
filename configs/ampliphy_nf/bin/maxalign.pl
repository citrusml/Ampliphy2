#!/usr/bin/perl
# Improves alignment by removing disturbing sequences.
# Pure perl version - no modules required.
#
# The basic idea is that in order to improve alignment some sequences
# can be removed from the input. "alignment area" is calculated as
# the number of free columns (columns with no gaps) multiplied with
# the number of sequences. To improve the score by removing sequences
# it is clear that the number of free columns must increase. Also, it
# is obvious that in order to get a new free column you must remove
# all sequences that has a gap in this column.
# NOTICE: Removing sets is which sequences to remove in order to
# improve the alignment area. Eliminating sets is which set of
# sequences it does not pay of to consider when removing sets.
#
# Copyright (C) 2007  Peter Wad Sackett, pws@cbs.dtu.dk
#                     Rodrigo de Oliveira, rodrigo@cbs.dtu.dk
# This program is free software under the GNU General Public License
# See "maxalign.pl -V" for details
#
use strict;

# CONSTANTS
my $reportfile = 'report.txt';
my $heuristicfastafile = 'heuristic.fsa';
my $heuristicincludeheaders = 'heuristic_include_headers.txt';
my $heuristicexcludeheaders = 'heuristic_exclude_headers.txt';
my $optimalfastafile = 'optimal.fsa';
my $optimalincludeheaders = 'optimal_include_headers.txt';
my $optimalexcludeheaders = 'optimal_exclude_headers.txt';

# GLOBAL VARIABLES

# The original sequences that must have their alignment area 
# improved are stored here with their headers
my @sequenceheaders;
my @sequences;
# The name of the fasta file with the sequences
my $inputfastafile = '';
# The name of the second shadow fasta file with the sequences
my $shadowfastafile = '';
my @shadowsequences;
my @shadowheaders;
# How the sequences in the two files correspond
my @shadowlinks;
# The length of the longest sequence, smaller sequences will be elongated
# with gaps. Therefore it is best that all sequences has the same length.
my $longestseq = 0;
# The original alignment area and the working alignment area.
my $orig_alignmentarea = -1;
my $alignmentarea = 0;
my $heuristic_alignmentarea = 0;
# The gap matrix constructed from sequences (and working matrix). Basically
# 0's are put in filled positions and 1's in the gaps.
my @orig_gapmatrix;
my @gapmatrix;
# The current number of free columns
my $freecolumns;
my $orig_freecolumns = -1;
# The current number of columns with gaps that WILL be in the final result.
my $gapcolumns;
# The original number of sequences and a working copy
my $orig_seqno;
my $seqno;
# Sequences which MUST be in the final output, even though it will be non-optimal.
my %mustkeepseq;
# The union gap pattern of kept sequences
my $keepgappattern = '';
# A time measurement variable.
my $timethen = time;
# The sequence sets that has a gap in this column. The original is a string
# of 1's and 0's while the other is a true bit vector for performance.
my @orig_sets;
my @sets;
# Where the gaps in the sets are.
my @gaps;
# The number of iterations the heuristic should go through.
my $iterations = -1;
# The number it had to go through in order to get highest alignment area.
my $iterationscount = 0;
# The standard debug level
my $debug = 0;
# Web output
my $www = 0;
# The heuristically excluded sequences
my %excludeseq;
# Translation table between orig_gapmatrix and working gapmatrix (and sets).
# Basically All-seqs - %excludeseqs
my @translation;
# Heuristic datapoints
my @heu_ex;
my @heu_alignvol;
# The optimal solution(s if equal)
my @solutions;
# The method the heuristic uses for deciding which sets to remove
# See sub: greatest_impact_set
my $heuristic_method = 1;
# Prove optimality - run branch and bound
my $optimality = 0;
# Detailed report
my $details = 0;
# Prefer piping
my $piping = 0;
# Default timeout in seconds for branch and bound
my $timeout = 300;
# Remove columns with all gaps in output
my $allgapsremoval = 0;
#### MAIN PROGRAM ####

&commandline_parsing;
&reading_fastafile;
&create_matrices;
   
&benchmark('Initialisation');

# The heuristic is of the greedy type; find the set of sequences which
# removal improves the alignment area most - remove the sequences.
# Rinse and repeat until no improvement can be made.
#### HEURISTIC ITERATION LOOP ####
while ($iterations > $iterationscount or $iterations == -1) {
   $iterationscount++;
   print "Iteration: $iterationscount\n" if $debug > 1;
   # Creating sets for this iteration
   &set_creation;
   # Joining/eliminating congruent sets and trivial too large sets.
   &congruent_set_joining;
   # Subset joining.
   &subset_joining;
   # Find best set to eliminate
   my ($bestset, $newalignvol) = &greatest_impact_set;
   if ($alignmentarea >= $newalignvol) {
      # Done
      $iterationscount--;
      print "No possible improvement on alignment area\n\n" if $debug > 1;      
      last; }
   print "New alignment area: $newalignvol\n\n" if $debug > 1;
   # Ready next round of heuristic
   $alignmentarea = $newalignvol;
   my $bitstring = unpack('b*',$sets[$bestset]);
   my @exseq;
   my $pointer = 0;
   while (-1 != ($pointer = index($bitstring, '1', $pointer))) {
      $excludeseq{$translation[$pointer]} = 1;
      push(@exseq, $translation[$pointer]);
      $pointer++; }
   @translation = ();
   for (my $i = 0; $i <= $#sequences; $i ++) {
      push(@translation, $i) unless exists $excludeseq{$i}; }
   # Store heuristic datapoints
   push(@heu_alignvol, $newalignvol);
   push(@heu_ex, [@exseq]);
} # End of heuristic loop

&optiontransformation(0);
unless ($piping) {
   &output_heuristic_data; }
elsif (not $optimality) {
   &output_piping;
   exit; }

print "Heuristic alignment area: $alignmentarea\n" if $debug;
&benchmark('Heuristic');

exit unless $optimality; 

# Re-establish original sequences and sets
print "Preparing for optimality proof\n" if $debug;
@translation = ();
for (my $i = 0; $i <= $#sequences; $i++) {
   push(@translation, $i); }
$heuristic_alignmentarea = $alignmentarea;
# Creating sets for Branch and Bound
&set_creation;
# Joining/eliminating congruent sets and trivial too large sets.
&congruent_set_joining;
# Subset joining.
&subset_joining;

&set_elimination;

&output_init_optimize unless $piping;

#### BEGINNING OF BRANCH AND BOUND ####
# Testing all set combinations.
# All set combinations should be tried now in order to find the combination
# that gives the highest alignment score. This can quickly be a very large
# time consuming task. Therefore the more sets we can eliminate before this
# step the better.
# It can be an immense performance improvement of the branch and bound
# algorithm if the sets can be ordered in some fashion that will force
# early pruning of the search tree. 4x has been observed.
# FAILED: Ordering sets by impact on the alignment area.
# FAILED: Ordering sets by their size.
# SUCCESS: Ordering sets by how many other sets they don't like.
# Sets can dislike (i.e. unable to both be chosen for removal) for a
# number of reasons; They can be subsets/supersets of the other, if
# both are removed they constitute a "too large set" definition.
my @dislikes;
my $lastdislikeset = -1;
&set_reordering;

# A few bits of optimizing
my $setsno = scalar @sets;
#my $maxremoveset = $iterationscount + 2; # used only for one failed bound algorithm
my $complete = 0;  # How many end nodes reached in search tree
my $pruning = 0;   # How many prunings
# A problem is an array of (problemstring, score, pointer, unionsets, unionsgaps)
my @stack = (['X' x $setsno, 0, pack('b*', '0' x $seqno), pack('b*', '0' x $longestseq)]);
# Set the alarm to stop after a given time
my $timeisup = 0;
$SIG{'ALRM'} = sub {$timeisup = 1; $SIG{'ALRM'} = 'IGNORE';};
alarm($timeout) if $timeout;
use integer; # Will save a few secs
while (scalar @stack) {
   last if $timeisup;
   # Get next problem in stack
   my ($problem, $pointer, $unionsets, $uniongaps) = @{pop(@stack)};
   for (;;) {
      # Bounding algorithm
      # It is vitally important for a successfull branch and bound algorithm
      # that we are able to calculate a good upper bound, i.e. a prediction
      # that lies within 2% of the best real alignment area possible for the
      # subproblem we are considering, but NEVER less than the best possible
      # value for the subproblem.
      # This calculation is simple but EXTREMELY efficient in pruning - remove to test.
      # The idea is that you calculate the positive effects of removing all the rest
      # of the sets in the subproblem, i.e. how many columns will be become ungapped,
      # while not including the negative effects, i.e. less sequences in the alignment.
      # The you can calculate the ideal (best possible) alignment area based on
      # the decisions already made and if that is less than what you know you can
      # achieve, then prune.
      my $testuniongaps = $uniongaps;
      my $j = $pointer;
      while (-1 != ($j = index($problem, 'X', $j))) {
          $testuniongaps |= $gaps[$j];
	  $j++; }
      if (($freecolumns + unpack('%32b*', $testuniongaps)) * ($seqno - unpack('%32b*', $unionsets)) < $alignmentarea) {
        $pruning++;
	last; }
      # Find next valid set to branch on among the constraints built in
      # to improve early pruning.
      while (-1 != ($pointer = index($problem, 'X', $pointer))) {
         # If this set is a subset of the union of already chosen sets
         # then it is automatically included (with gaps) for performance.
         last if unpack('%32b*', ($unionsets | $sets[$pointer]) ^ $unionsets);
         $uniongaps |= $gaps[$pointer];
         substr($problem, $pointer, 1, '1');
         $pruning++;
         $pointer++; }

      if ($pointer != -1) {
         # Branching on a subproblem
         substr($problem, $pointer, 1, '0');
         push(@stack, [$problem, $pointer+1, $unionsets, $uniongaps]);
         substr($problem, $pointer, 1, '1');
         $unionsets |= $sets[$pointer];
         $uniongaps |= $gaps[$pointer];
	 # Pruning - cuts running time with 25%
         my $setsize = unpack('%32b*', $unionsets);
	 my $gapno = unpack('%32b*', $uniongaps);
         if ($alignmentarea > ($longestseq - $gapcolumns) * ($seqno - $setsize) or
              ($freecolumns*$setsize) > $gapno*($seqno-$setsize)) {
            $pruning++;
	    last; }
	 # OK problem - back to work
         if ($pointer < $lastdislikeset) {
            foreach my $badsets (@{$dislikes[$pointer]}) {
               substr($problem, $badsets, 1, '0'); } }
         $pointer++; 
         next; }

      # We have reached an end node of the search tree. Is it a solution?
      $complete++;
      my $score = ($freecolumns + unpack('%32b*', $uniongaps)) * ($seqno - unpack('%32b*', $unionsets));
      last if $score < $alignmentarea; # No
      if ($score == $alignmentarea) {
	 print "Equal solution found: $problem\n" if $debug > 1; 
         push(@solutions, $unionsets);
         last; }
      print "Better solution found: $problem - $score\n" if $debug > 1;
      @solutions = ($unionsets);
      $alignmentarea = $score;
      # IMPORTANT IDEA: Possible eliminations of sets here,
      # but only if the heuristic is bad, which it isn't
      last;
   }
} # End of branch and bound

unless ($timeisup) {
   if ($heuristic_alignmentarea != $alignmentarea) {
      # Choosing solution with fewest removed sequences
      @solutions = sort {unpack('%32b*', $a) <=> unpack('%32b*', $b)} @solutions if $#solutions > 0; 
      %excludeseq = ();
      my $bitstring = unpack('b*',$solutions[0]);
      my $pointer = 0;
      while (-1 != ($pointer = index($bitstring, '1', $pointer))) {
         $excludeseq{$pointer} = 1;
         $pointer++; } 
      &optiontransformation(1);  } }

if ($debug) {
   if ($timeisup) {
      print "Program ran out of time trying to prove optimality.\n"; }
   else {
      print "Optimal alignment area: $alignmentarea\n"; } }

if ($piping) {
   &output_piping; }
else {
   &output_end_optimize; }

&benchmark('Branch and Bound');


#### SUBROUTINES ####


sub usage {
   my ($msg) = @_;
   print "$msg\n\n" if defined $msg;
   print "Usage: maxalign.pl [-a] [-d] [-f?] [-h?] [-i?] [-o[?]] [-v?] [-w] [<fastafile>] [<fastafile>]\n";
   print "   -a  All gaps columns in output is removed\n";
   print "   -d  Detailed report\n";
   print "   -f? Prepends ? to output files. This could be a path and/or filename addition.\n";
   print "   -h? Heuristic to use. Replace ? with 1-3. Default is 1 and that seems optimal.\n";
   print "   -i? Iteratively removes up til ? sets until no improvement (0-999)\n";
   print "   -o? Prove optimality in max ? minuts. 0 means no limit. Default is 5 min.\n";
   print "   -p  Pipe only the resulting alignment to STDOUT. No files - no fuss.\n";
   print "   -v? Verbosity level (1-3)\n";
   print "   -w  Web output\n";
   print "   -V  Version and copyright\n";
   print "Options must be given separately. If a fasta file is not given on command line,\n";
   print "then the alignment is expected as fasta on STDIN.\n";
   print "IFF two fasta files are given on command line then the optimal alignment area is\n";
   print "calculated for the first file, and the disturbing sequences from this calculation will\n";
   print "be removed from the second file, which would be output, i.e. the calculation will be\n";
   print "based on the first file, the effect on the second. This requires that the fasta IDs\n";
   print "are identical in both files. Used with amino acid and nucleotide fasta files - one each.\n";
   print "If a plus, +, is the first character in the sequence header (>+AC0023),\n";
   print "then that sequence is always included in the output. The + will be removed.\n";
   print "Possibly produces the output files; $reportfile, $heuristicfastafile, $heuristicincludeheaders\n";
   print "$heuristicexcludeheaders, $optimalfastafile, $optimalincludeheaders and $optimalexcludeheaders\n";
   exit;
}


sub commandline_parsing {
   # Parsing options
   while (scalar @ARGV) {
      if ($ARGV[0] =~ m/^-a$/) {
         $allgapsremoval = 1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-d$/) {
         $details = 1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-f(.+)$/) {
         $reportfile = $1 . $reportfile;
         $heuristicfastafile = $1 . $heuristicfastafile;
         $heuristicincludeheaders = $1 . $heuristicincludeheaders;
         $heuristicexcludeheaders = $1 . $heuristicexcludeheaders;
         $optimalfastafile = $1 . $optimalfastafile;
         $optimalincludeheaders = $1 . $optimalincludeheaders;
         $optimalexcludeheaders = $1 . $optimalexcludeheaders;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-h([123])$/) {
         $heuristic_method = $1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-i(\d{1,3})$/) {
         $iterations = $1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-o(\d*)$/) {
         $optimality = 1;
	 $timeout = 60 * $1 if defined $1 and $1 ne '';
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-p$/) {
         $piping = 1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-v([123])$/) {
         $debug = $1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-w$/) {
         $www = 1;
         shift @ARGV; }
      elsif ($ARGV[0] =~ m/^-V$/) {
         &version_copyright; }
      elsif ($ARGV[0] =~ m/^-/) {
         &usage('Unknown option'); }
      elsif ($inputfastafile eq '') {
         $inputfastafile = shift @ARGV; }
      elsif ($shadowfastafile eq '') {
         $shadowfastafile = shift @ARGV; }
      else {
         &usage(); } }
   &usage("Option -p can't be used together with -d, -v or -w") if $piping and ($details or $debug or $www);
}


# Reading input file, sequence headers go to @sequenceheaders
# The corresponding sequences go to @sequences. The number of
# sequences ($seqno) and the length of the longest sequence 
# ($longestseq) is also calculated.
# Reads and verifies the extra shadow fastafile.
sub reading_fastafile {
   # Normal fasta file
   my $IN;
   if ($inputfastafile) {
      open($IN, $inputfastafile) or &usage("Can't read file $inputfastafile; reason $!"); }
   else {
      $IN = *STDIN; }
   my ($header, $seq, $lineno, $line) = ('', '', 0);
   eval {
      local $SIG{'ALRM'} = sub { die "alarm\n" };
      alarm(1);
      $line = <$IN>;
      alarm(0); };
   if ($@) {
      die "$@\n" if $@ ne "alarm\n";
      &usage('Input on STDIN is expected, terminating MaxAlign.'); }  
   while (defined $line ) {
      $lineno++;
      chomp $line;
      if ($line =~ m/^>(.+)/) {
         if ($header) {
  	    if (substr($header, 0, 1) eq '+') {
	       substr($header, 0, 1, '');
	       my $no = scalar @sequenceheaders;
	       $mustkeepseq{$no} = 1; }
            push(@sequenceheaders, $header);
   	    $seq =~ s/\s//g;
	    $longestseq = length($seq) if $longestseq < length($seq);
	    push(@sequences, $seq); }
         $header = $1;
         $seq = ''; }
      elsif ($line =~ m/[^\w\s\-]/) {
         &report_user_error($line) unless $piping;
         &usage("Bad line in input line $lineno:\n $line"); }
      else {
         $seq .= $line; }
      $line = <$IN>; }
   close $IN;
   if ($header) {
      # FIX: check/strip '+' BEFORE push to correctly identify must-keep sequence and prevent trailing '+'
      if (substr($header, 0, 1) eq '+') {
	 substr($header, 0, 1, '');
	 my $no = scalar @sequenceheaders;
	 $mustkeepseq{$no} = 1; }
      push(@sequenceheaders, $header);
      $seq =~ s/\s//g;
      $longestseq = length($seq) if $longestseq < length($seq);
      push(@sequences, $seq); }
   $orig_seqno = scalar @sequenceheaders;
   unless ($orig_seqno) {
      &report_user_error() unless $piping;
      &usage("Input is not fasta format"); }
   # Shadow fasta file read
   if ($shadowfastafile) {
      my %ids;
      for (my $i = 0; $i <= $#sequenceheaders; $i++) {
         $sequenceheaders[$i] =~ m/^(\S+)/;
	 if (exists $ids{$1}) {
            &report_user_error("$i $sequenceheaders[$i] $1 - ID is not unique") unless $piping;
            &usage("$1 - ID is not unique (first file):\n"); }
	 $ids{$1} = 0;
	 $shadowlinks[$i] = $1; }
      open($IN, $shadowfastafile) or &usage("Can't read file $shadowfastafile; reason $!");
      ($header, $seq, $lineno, $line) = ('', '', 0);
      $line = <$IN>;
      while (defined $line ) {
         $lineno++;
         chomp $line;
         if ($line =~ m/^>(.+)/) {
            if ($header) {
   	       substr($header, 0, 1, '') if substr($header, 0, 1) eq '+';
               push(@shadowheaders, $header);
   	       $seq =~ s/\s//g;
	       push(@shadowsequences, $seq); }
            $header = $1;
	    $header =~ m/^(\S+)/;
	    unless (exists $ids{$1}) {
               &report_user_error("$line\nID $1 does not exists in other file") unless $piping;
               &usage("ID $1 does not exists in other file:\n $line"); }
	    elsif ($ids{$1}) {
               &report_user_error("$line\nID $1 is not unique") unless $piping;
               &usage("ID $1 is not unique:\n $line"); }
	    else {
	       $ids{$1} = 1; }
            $seq = ''; }
         elsif ($line =~ m/[^\w\s\-]/) {
            &report_user_error($line) unless $piping;
            &usage("Bad line in input line $lineno:\n $line"); }
         else {
            $seq .= $line; }
         $line = <$IN>; }
      close $IN;
      if ($header) {
         substr($header, 0, 1, '') if substr($header, 0, 1) eq '+';
         push(@shadowheaders, $header);
         $seq =~ s/\s//g;
         push(@shadowsequences, $seq); }
      if ($#sequences != $#shadowsequences) {
         &report_user_error("There is not the same amount of sequences in both files") unless $piping;
         &usage("There is not the same amount of sequences in both files"); }
      # Correspondence bewteen sequences in the files, for later use in output.
      for (my $i = 0; $i <= $#shadowheaders; $i++) {
         $shadowheaders[$i] =~ m/^(\S+)/;
	 $ids{$1} = $i; }
      for (my $i = 0; $i <= $#shadowlinks; $i++) {
         $shadowlinks[$i] = $ids{$shadowlinks[$i]}; }
      my $longseq = 0;
      for (my $i = 0; $i <= $#shadowsequences; $i++) {
	 $longseq = length $shadowsequences[$i] if $longseq < length $shadowsequences[$i]; }
      for (my $i = 0; $i <= $#shadowsequences; $i++) {
	 $shadowsequences[$i] .= '-' x  ($longseq - length $shadowsequences[$i])
	     if $longseq > length $shadowsequences[$i]; } 
   }
}


# The gap matrix is created. All sequences are transformed into
# a matrix where a row is a sequence where amino acids are replaced
# with a 0 and a gap is replaced with 1. If some sequences are
# "short" then they are filled with gaps.
# After that a set matrix is created.
sub create_matrices {
   # Gap matrix
   for (my $i = $#sequences; $i >= 0; $i--) {
      my $seq = $sequences[$i];
      $seq .= '-' x ($longestseq - length($seq)) if $longestseq > length($seq);
      $seq =~ tr/\-/0/c;
      $seq =~ tr/\-/1/;
      $orig_gapmatrix[$i] = [split(//, $seq)]; }
   # Sets (matrix)
   for (my $i = 0; $i < $longestseq; $i++) {
      my $bitstring = '';
      for (my $j = 0; $j <= $#orig_gapmatrix; $j++) {
         $bitstring .= $orig_gapmatrix[$j][$i]; }
      push(@orig_sets, $bitstring); }
   # Create the gap pattern that is certain from kept sequences
   $keepgappattern = '0' x $longestseq;
   foreach my $i (keys %mustkeepseq) {
      for (my $j = 0; $j < $longestseq; $j++) {
         substr($keepgappattern, $j, 1, '1') if $orig_gapmatrix[$i][$j] eq '1'; } }
   # Initialising translation table
   for (my $i = 0; $i <= $#sequences; $i++) {
      push(@translation, $i); }
}


# Creating sets. A set (@sets) must be created for each amino acid
# position of the current working sequences that has a gap on that position.
# These sets are subsets of the original input. If there are no gaps
# in any sequence at that position, then it is ignored, since no
# improvement can be made. It also creates a gap vector (@gaps)
# for each set that denotes the column in question. The alignment
# area is also calculated for the input.
sub set_creation {
   # Preparing new gap matrix, can actually be skipped, since it is
   # not used in the current code, but it is quick, so ...
   #@gapmatrix = ();
   #for (my $i = $#translation; $i >= 0; $i--) {
   #   $gapmatrix[$i] = $orig_gapmatrix[$translation[$i]]; }
   #if ($debug > 2) {
   #   print "Gap matrix\n";
   #   for (my $i = 0; $i <= $#gapmatrix; $i++) {
   #      for (my $j = 0; $j < $longestseq; $j++) {
   #         print $gapmatrix[$i][$j]; }
   #      print "\n"; } }
   # Set creation
   @sets = ();
   @gaps = ();
   my @exseq = ();
   @exseq = keys %excludeseq if scalar @translation != scalar @sequences;
   for (my $i = 0; $i < $longestseq; $i++) {
      next if substr($keepgappattern, $i, 1) eq '1';
      my $bitstring = $orig_sets[$i];
      foreach my $sq (@exseq) {
         substr($bitstring, $sq, 1, 'X'); }
      next unless $bitstring =~ tr/1/1/;
      $bitstring =~ tr/01X/01/d;
      push(@sets, pack('b*', $bitstring));
      $bitstring = '0' x $longestseq;
      substr($bitstring, $i, 1, '1');
      push(@gaps,pack('b*', $bitstring)); }
   # New values
   $seqno = scalar @translation;
   $freecolumns = $longestseq - scalar @sets + ($keepgappattern =~ tr/1/1/);
   $orig_freecolumns = $freecolumns if $orig_freecolumns < 0;
   $alignmentarea = $freecolumns * $seqno if $alignmentarea < $freecolumns * $seqno;
   $orig_alignmentarea = $alignmentarea if $orig_alignmentarea < 0;
   # Debug code
   if ($debug > 1) {
      print "Number of sequences: $seqno\n";
      print "Longest sequence: $longestseq\n";
      print "Free columns: $freecolumns\n";
      print "Original alignment area: $orig_alignmentarea\n"; 
      print "Sets after creation: ", scalar @sets, "\n"; }
   if ($debug > 2) {
      print "Sets\n";
      for (my $i = 0; $i <= $#sets; $i++) {
         my $seq = unpack('b*',$sets[$i]);
         print substr($seq, 0, $seqno), "\n"; } }
}


# Joining/eliminating congruent sets and too large sets.
# In order to maximize the alignment area, one must figure out
# which sets to remove - a difficult task.
# Some sets are identical to each other, i.e. two sequences that
# has two gaps in the same position. These sets and their gap 
# vectors are combined thereby eliminating a set, making the task
# of selecting sets easier.
# Some sets are so large (contain so many sequences) that if they
# are removed, that no matter how well we free up other columns
# we can never improve the original alignment area. These sets
# also eliminated from further consideration.
sub congruent_set_joining {
   # First find the set size order
   my %setorder;
   my @setsize;
   for (my $i = $#sets; $i >= 0; $i--) {
      $setorder{$i} = $setsize[$i] = unpack('%32b*', $sets[$i]); }
   my @order = sort {$setorder{$a} <=> $setorder{$b}} keys %setorder;
   # Sets are congruent if they are the same size (performance) and
   # the same sequences (necessary condition).
   $gapcolumns = 0;
   for (my $i = 0; $i <= $#order; $i++) {
      my $orderi = $order[$i];
      if ($alignmentarea > $longestseq * ($seqno - $setsize[$orderi])) {
         $setsize[$orderi] = 0;  # Large set
         $gapcolumns++;
         next; }
      for (my $j = $i+1; $j <= $#order; $j++) {
         last if $setsize[$orderi] != $setsize[$order[$j]];
         next if unpack('%32b*', $sets[$orderi] ^ $sets[$order[$j]]);
         $gaps[$order[$j]] |= $gaps[$orderi];
         $setsize[$orderi] = 0; # Congruent set
         last; } }
   # Removing congruent eliminated and too large sets
   for (my $i = $#setsize; $i >= 0; $i--) {
      next if $setsize[$i];
      splice(@gaps, $i, 1);
      splice(@sets, $i, 1); }
   # Debugging
   if ($debug > 1) {
      print "Sets after eliminating congruent and large sets: ", scalar @sets, "\n";
      print "Minimum gapped columns: $gapcolumns\n"; }
   if ($debug > 2) {
      print "Sets and corresponding gap vectors\n";
      for (my $i = 0; $i <= $#sets; $i++) {
         my $seq = unpack('b*',$sets[$i]);
         print "$i, ", substr($seq, 0, $seqno), "\n";
         $seq = unpack('b*',$gaps[$i]);
         print "$i, ", substr($seq, 0, $longestseq), "\n"; } }
}


# Subset joining - adding sets to larger sets containing them.
# Some sets are subsets of other sets. So if we remove a larger set
# we may also have removed a smaller set and thereby freed more
# columns. Here the gap vectors of subsets are added to the larger set
# containing them, so this can be calculated.
sub subset_joining {
   # First find the set size order
   my %setorder;
   for (my $i = $#sets; $i >= 0; $i--) {
      $setorder{$i} = unpack('%32b*', $sets[$i]); }
   my @order = sort {$setorder{$a} <=> $setorder{$b}} keys %setorder;
   # Trivial joining if subset
   for (my $i = $#order; $i > 0; $i--) {
      my $orderi = $order[$i];
      for (my $j = $i-1; $j >= 0; $j--) {
         $gaps[$orderi] |= $gaps[$order[$j]]
            unless unpack('%32b*',($sets[$orderi] | $sets[$order[$j]]) ^ $sets[$orderi]); } }
   # Debugging
   if ($debug > 2) {
      print "Sets and corresponding gap vectors after joining\n";
      for (my $i = 0; $i <= $#sets; $i++) {
         my $seq = unpack('b*',$sets[$i]);
         print "$i, ", substr($seq, 0, $seqno), "\n";
         $seq = unpack('b*',$gaps[$i]);
         print "$i, ", substr($seq, 0, $longestseq), "\n"; } }
}


# Find and return the set (and impact) with the greatest the impact on the alignment area.
# How is impact measured? Here is an incomplete list of possibilities in descending effectiveness:
# 1) Most bang for the buck - best improvement per sequence removed
# 2) Greatest improvement on alignment area
# 3) Fewest sequences removed while still improving
sub greatest_impact_set {
   my $set = 0;
   my $impact = -1;
   if ($heuristic_method == 1) {
      my $bang = -1;
      for (my $i = $#sets; $i >= 0; $i--) {
         my $thisimpact = ($seqno - unpack('%32b*', $sets[$i])) * ($freecolumns + unpack('%32b*', $gaps[$i]));
         my $thisbang = ($thisimpact - $alignmentarea)/unpack('%32b*', $sets[$i]);
         next if $thisbang < $bang;
         if ($thisbang > $bang) {
            $bang = $thisbang;
            $impact = $thisimpact;
            $set = $i;
            next; }
         next if unpack('%32b*', $gaps[$i]) >= unpack('%32b*', $gaps[$set]);
         $impact = $thisimpact;
         $set = $i; } }
   elsif ($heuristic_method == 2) {
      for (my $i = $#sets; $i >= 0; $i--) {
         my $thisimpact = ($seqno - unpack('%32b*', $sets[$i])) * ($freecolumns + unpack('%32b*', $gaps[$i]));
         next if $thisimpact < $impact;
         if ($thisimpact > $impact) {
            $impact = $thisimpact;
            $set = $i;
            next; }
         next if unpack('%32b*', $gaps[$i]) >= unpack('%32b*', $gaps[$set]);
         $impact = $thisimpact;
         $set = $i; } }
   elsif ($heuristic_method == 3) {
      my $fewest = -1;
      for (my $i = $#sets; $i >= 0; $i--) {
         my $thisimpact = ($seqno - unpack('%32b*', $sets[$i])) * ($freecolumns + unpack('%32b*', $gaps[$i]));
         my $thisfewest = unpack('%32b*', $sets[$i]);
         next if $thisimpact <= $alignmentarea;
         next if $thisfewest < $fewest;
         if ($thisfewest > $fewest) {
            $fewest = $thisfewest;
            $impact = $thisimpact;
            $set = $i;
            next; }
         next if unpack('%32b*', $gaps[$i]) >= unpack('%32b*', $gaps[$set]);
         $impact = $thisimpact;
         $set = $i; } }
   else {
      die "Unknown heuristic method"; }
   return ($set, $impact);
}


# Now we can eliminate sets that are important to keep or too large.
# For every set eliminated this way we have
# actually calculated that we HAVE to have at least 1 gap column.
sub set_elimination {
   for (my $i = $#sets; $i >= 0; $i--) {
      my $setsize = unpack('%32b*', $sets[$i]);
# THIS IS NO GOOD
#      if (($freecolumns*$setsize) > unpack('%32b*', $gaps[$i])*($seqno-$setsize)) {
#         # Some sets are too important to remove; If the minimum gain of having
#         # them included, is greater than the maximum gain of removing them, then eliminate.
#         splice(@sets, $i, 1);
#         splice(@gaps, $i, 1); }
#      els
      if ($alignmentarea > ($longestseq - $gapcolumns) * ($seqno - $setsize)) {
         # If a set is so large, that its removal under no circumstances
         # can improve alignment area, it can be eliminated.
         splice(@sets, $i, 1);
         splice(@gaps, $i, 1); } }
   $gapcolumns = &getgapcolumns;
   my $lastgapcolumns = -1;
   # Large size check is actually iterative
   while ($lastgapcolumns != $gapcolumns) {
      $lastgapcolumns = $gapcolumns;
      for (my $i = $#sets; $i >= 0; $i--) {
         next if $alignmentarea <= ($longestseq - $gapcolumns) * ($seqno - unpack('%32b*', $sets[$i]));
         splice(@gaps, $i, 1);
         splice(@sets, $i, 1); 
         $gapcolumns = &getgapcolumns; } }
   # Debugging
   if ($debug > 1) {
      print "Sets after eliminating too important sets: ", scalar @sets, "\n";
      print "New minimum gapped columns: $gapcolumns\n"; }
   if ($debug>2) {
      print "Sets and corresponding gap vectors\n";
      for (my $i = 0; $i <= $#sets; $i++) {
         my $seq = unpack('b*',$sets[$i]);
         print "$i, ", substr($seq, 0, $seqno), "\n";
         $seq = unpack('b*',$gaps[$i]);
         print "$i, ", substr($seq, 0, $longestseq), "\n"; } }
}


# By reordering the sets according to their dislikes and using that
# fact in implicit pruning considerably performace increase is achieved.
sub set_reordering {
   # First find out which sets dislikes each other
   my %setorder = ();
   for (my $i = $#sets; $i >= 0; $i--) {
      $setorder{$i} = {()}; }
   for (my $i = 0; $i < $#sets; $i++) {
      my $set_i = $sets[$i];
      for (my $j = $i+1; $j <= $#sets; $j++) {
         my $set_j = $sets[$j];
         # Dislike if one set is subset of the other
         my $union = unpack('b*', $set_i | $set_j);
         if ($union eq unpack('b*', $set_i) or $union eq unpack('b*', $set_j)) {
            $setorder{$i}{$j} = 1;
            $setorder{$j}{$i} = 1;
            next; }
	 my $setssize = unpack('%32b*', $set_i | $set_j);
         # Dislike if sets united is a "too large set"
         if ($alignmentarea > ($longestseq - $gapcolumns) * ($seqno - $setssize)) {
            $setorder{$i}{$j} = 1;
            $setorder{$j}{$i} = 1;
            next; }
	 # Dislike if sets united are "too important to remove"
# No good
#         if (($freecolumns*$setssize) > unpack('%32b*', $gaps[$i] | $gaps[$j])*($seqno-$setssize)) {
#            $setorder{$i}{$j} = 1;
#            $setorder{$j}{$i} = 1;
#            next; }
      } }
   # The reordering. Primary and secondary sorting is important.
   # Primary sorting; Number of dislikes - most is best
   # Secondary sorting; Set size - largest is best
   # Tertiary sorting; Number of gaps - most is best
   my @order = ();
   my %backward = ();
   my @backlist = ();
   for (my $i = 0; $i <= $#sets; $i++) {
      my @k = sort { my $res = scalar keys %{$setorder{$b}} <=> scalar keys %{$setorder{$a}};
                     return $res if $res;
                     $res = unpack('%32b*', $sets[$b]) <=> unpack('%32b*', $sets[$a]);
                     return $res if $res;
                     return unpack('%32b*', $gaps[$b]) <=> unpack('%32b*', $gaps[$a]);
                   } keys %setorder;
      my $winner = $k[0];
      push(@order, $winner);
      $backward{$winner} = $#order;
      foreach my $key (keys %{$setorder{$winner}}) {
         delete $setorder{$key}{$winner}; }
      push(@backlist, [keys %{$setorder{$winner}}]);
      delete $setorder{$winner}; }
   my @tmpgaps;
   my @tmpsets;
   for (my $i = 0; $i <= $#order; $i++) {
      push(@tmpgaps, $gaps[$order[$i]]);
      push(@tmpsets, $sets[$order[$i]]); }
   @gaps = @tmpgaps;
   @sets = @tmpsets;
   # Build up the dislikes list;
   for (my $i = 0; $i <= $#sets; $i++) {
      my @bad = ();
      foreach my $set_no (@{$backlist[$i]}) {
         push(@bad, $backward{$set_no}); }
      $lastdislikeset = $i if $lastdislikeset == -1 and scalar @bad == 0;
      @bad = sort {$a <=> $b} @bad if scalar @bad;
      push(@dislikes, [@bad]); }
   # The debugging bit
   if ($debug > 1) {
      for (my $i = 0; $i <= $#sets; $i++) {
         print "$i => @{$dislikes[$i]}\n"; } }
}


sub output_heuristic_data {
   # Included headers
   open(OUT, '>', $heuristicincludeheaders) or
      die "Can't write file: $heuristicincludeheaders\nReason: $!\n";
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT "$sequenceheaders[$i]\n" unless exists $excludeseq{$i}; }
   close OUT;
   # Best alignment fasta file for heuristic
   open(OUT, '>', $heuristicfastafile) or
      die "Can't write file: $heuristicfastafile\nReason: $!\n";
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT &fasta($i) unless exists $excludeseq{$i}; }
   close OUT;
   # Excluded headers
   open(OUT, '>', $heuristicexcludeheaders) or
      die "Can't write file: $heuristicexcludeheaders\nReason: $!\n";
   for (my $i = 0; $i <= $#heu_ex; $i++) {
      if ($details) {
         print OUT '# Iteration=', $i+1, ' alignmentarea=', $heu_alignvol[$i];
         print OUT ' ExcludedSeqs=', scalar @{$heu_ex[$i]}, "\n"; }
      foreach my $no (@{$heu_ex[$i]}) {
         print OUT "$sequenceheaders[$no]\n"; } }
   close OUT;
   # The report
   open(OUT, '>', $reportfile) or
      die "Can't write file: $reportfile\nReason: $!\n";
   print OUT &tag('Results from CBS MaxAlign Service', 'h3'), "\n\n";
   print OUT 'Original number of sequences: ', &tag($orig_seqno, 'b', 'br');
   print OUT 'Original total number of columns: ', &tag($longestseq, 'b', 'br');
   print OUT 'Original ungapped columns: ', &tag($orig_freecolumns, 'b', 'br');
   print OUT 'Original alignment area: ', &tag($orig_alignmentarea, 'b', 'br', 'br');
   if ($details) {
      print OUT &tag('Details from the heuristic', 'b', 'br');
      print OUT 'Method for sequence set removal: ';
      print OUT &tag('Bang for the buck', 'br', 'br') if $heuristic_method == 1;
      print OUT &tag('Greatest improvement', 'br', 'br') if $heuristic_method == 2;
      print OUT &tag('Fewest removed', 'br', 'br') if $heuristic_method == 3; }
   if ($iterationscount == 0) {
      print OUT &tag("This alignment can't be improved by removing sequences");
      close OUT;
      exit; }
   if ($details) {
      my $accuxno = 0;
      for (my $i = 0; $i <= $#heu_ex; $i++) {
         print OUT &tag('Iteration: ' . ($i+1), 'br');
         my $xno = scalar @{$heu_ex[$i]};
         $accuxno += $xno;
         print OUT &tag("Excluded sequences this round: $xno", 'br');
         my $freecols = $heu_alignvol[$i] / ($seqno-$accuxno);
         print OUT &tag("Ungapped columns: $freecols", 'br');
         print OUT &tag("New alignment area: $heu_alignvol[$i]", 'br', 'br'); } }
   if ($details) {
      print OUT 'Number of sequences in heuristic result: ', &tag($orig_seqno-scalar keys %excludeseq, 'b', 'br');
      print OUT 'Total columns in heuristic result: ', &tag(length($sequences[0]), 'b', 'br') if $allgapsremoval;
      print OUT 'Ungapped columns in heuristic result: ', &tag($freecolumns, 'b', 'br');
      print OUT 'Heuristic best alignment area: ', &tag($heu_alignvol[-1], 'b', 'br', 'br');
      print OUT 'Resulting included headers from heuristic: ', &tag($heuristicincludeheaders, 'a', 'br');
      print OUT 'Resulting excluded headers from heuristic: ', &tag($heuristicexcludeheaders, 'a', 'br', 'br'); }
   else {
      print OUT 'Number of sequences in result: ', &tag($orig_seqno-scalar keys %excludeseq, 'b', 'br');
      print OUT 'Total columns in heuristic result: ', &tag(length($sequences[0]), 'b', 'br') if $allgapsremoval;
      print OUT 'Ungapped columns in result: ', &tag($freecolumns, 'b', 'br');
      print OUT 'Resulting alignment area: ', &tag($heu_alignvol[-1], 'b', 'br', 'br');
      print OUT 'Resulting included headers: ', &tag($heuristicincludeheaders, 'a', 'br'); 
      print OUT 'Resulting excluded headers: ', &tag($heuristicexcludeheaders, 'a', 'br', 'br'); }
   print OUT 'Download your resulting alignment: ', &tag($heuristicfastafile, 'a', 'br');
   print OUT &tag('or copy-paste from below. ', 'br'); 
   &output_alignment unless $optimality;
   close OUT;
}


sub output_init_optimize {
   open(OUT, '>>', $reportfile) or
      die "Can't write file: $reportfile\nReason: $!\n";
   print OUT &tag('','br');
   print OUT &tag('Details from the proof-of-optimality algorithm', 'b', 'br');
   print OUT &tag('Number of sets: ' . scalar @sets, 'br');
   print OUT &tag('Theoretical number of solutions: 2^' . scalar @sets, 'br');
   close OUT;
}


sub output_end_optimize {
   # Report
   open(OUT, '>>', $reportfile) or
      die "Can't write file: $reportfile\nReason: $!\n";
   if ($timeisup) {
      print OUT &tag('The program was not given enough computer time to prove optimality.', 'b', 'br');
      print OUT &tag("Meanwhile, here is the result from the heuristic.", 'br');
      &output_alignment;
      close OUT;
      exit; }
   if (scalar @solutions < 1) {
      print OUT &tag('An error has occurred since no solutions are found', 'br');
      print OUT &tag('and it was possible to find a solution via the heuristic.', 'br');
      print OUT &tag('Please, contact the author and report the problem and your input.', 'br', 'br');
      print OUT &tag("Meanwhile, here is the result from the heuristic", 'br');
      &output_alignment;
      close OUT;
      exit; }
   if ($heuristic_alignmentarea == $alignmentarea) {
      print OUT &tag('This solution is proved to be optimal', 'b', 'br') ; }
   else {
      print OUT &tag('Above solution is NOT optimal', 'b', 'br') ; 
      print OUT &tag('Here is a better solution', 'br') ; }
   if ($details or  $heuristic_alignmentarea != $alignmentarea) {
      print OUT 'Number of optimal solutions: ', &tag(scalar @solutions, 'br');
      if (scalar @solutions > 1) {
         print OUT &tag('Choosing the solution with fewest sequences eliminated', 'br'); }
      print OUT &tag("Number of solutions checked: $complete", 'br') if $details;
      print OUT &tag("Number of explicit prunings in search tree: $pruning", 'br') if $details;
      print OUT 'Number of sequences in optimized result: ', &tag($seqno-unpack('%32b*', $solutions[0]), 'b', 'br');
      print OUT 'Total columns in optimized result: ', &tag(length($sequences[0]), 'b', 'br') if $allgapsremoval;
      print OUT 'Ungapped columns: ', &tag($alignmentarea/($seqno-unpack('%32b*', $solutions[0])), 'b', 'br');
      print OUT 'Optimal alignment area: ', &tag($alignmentarea, 'b', 'br');
      print OUT 'Resulting fasta file from optimizing: ', &tag($optimalfastafile, 'a', 'br');
      print OUT 'Resulting included headers from optimizing: ', &tag($optimalincludeheaders, 'a', 'br');
      print OUT 'Resulting excluded headers from optimizing: ', &tag($optimalexcludeheaders, 'a', 'br');
      print OUT &tag('', 'br'), &tag('Thanks for visiting - Come again soon', 'br') if $www; }
   &output_alignment;
   close OUT; 
   my $res = unpack('b*', $solutions[0]);
   # Optimal fasta file
   open(OUT, '>', $optimalfastafile) or
      die "Can't write file: $optimalfastafile\nReason: $!\n";
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT &fasta($i) if substr($res, $i, 1) eq '0'; }
   close OUT;
   # Optimal included headers
   open(OUT, '>', $optimalincludeheaders) or
      die "Can't write file: $optimalincludeheaders\nReason: $!\n";
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT "$sequenceheaders[$i]\n" if substr($res, $i, 1) eq '0'; }
   close OUT;
   # Optimal excluded headers
   open(OUT, '>', $optimalexcludeheaders) or
      die "Can't write file: $optimalexcludeheaders\nReason: $!\n";
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT "$sequenceheaders[$i]\n" if substr($res, $i, 1) eq '1'; }
   close OUT;
}


sub output_alignment {
   print OUT &tag('New alignment in fasta format:', 'br', 'br');
   print OUT "<pre>\n" if $www;
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print OUT &fasta($i) unless exists $excludeseq{$i}; }
   print OUT "</pre>\n" if $www;
}



sub output_piping {
   for (my $i = 0; $i < $orig_seqno; $i++) {
      print &fasta($i) unless exists $excludeseq{$i}; }
}


sub optiontransformation {
   my $level = shift @_;
   if ($shadowfastafile) {
      if ($level == 0) {
         @sequences = @shadowsequences;
         @sequenceheaders = @shadowheaders; }
      elsif ($allgapsremoval) {
         @sequences = @shadowsequences; }
      my %tmp;
      foreach my $key (keys %excludeseq) {
         $tmp{$key} = 1; }
      %excludeseq = %tmp; }
   if ($allgapsremoval) {
      my %positions;
      for (my $i = length($sequences[0]) - 1; $i >= 0; $i--) {
         $positions{$i} = 1; }
      for (my $i = 0; $i <= $#sequences; $i++) {
         next if exists $excludeseq{$i};
	 foreach my $keypos (keys %positions) {
	    delete $positions{$keypos} if substr($sequences[$i], $keypos, 1) ne '-'; } }
      foreach my $keypos (sort {$b <=> $a} keys %positions) {
         for (my $i = 0; $i <= $#sequences; $i++) {
            substr($sequences[$i], $keypos, 1, ''); } } }
}


sub report_user_error {
   my ($line) = @_;
   open(OUT, '>', $reportfile) or
      die "Can't write file: $reportfile\nReason: $!\n";
   print OUT &tag('Your data is not in fasta format', 'h3');
   print OUT &tag("The problem occurs in line $.:", 'br');
   print OUT &tag($line, 'br', 'br');
   print OUT 'You can see an example of fasta format at ',
                 &tag('http://www.cbs.dtu.dk/services/fasta.example.php', 'a', 'br');
}  


# Calculates the minimum number of gapped columns for the remaining sets
sub getgapcolumns {
   my $union = pack('b*', '0' x $longestseq);
   for (my $i = $#sets; $i >= 0; $i--) {
      $union |= $gaps[$i]; }
   return $longestseq - $freecolumns - unpack('%32b*', $union);
}
     

# Time measurement subroutine
sub benchmark {
   return unless $debug;
   my ($msg) = @_;
   my $timenow = time;
   print "$msg: ", $timenow - $timethen, " sec\n\n";
   $timethen = $timenow;
}


# Returns a nicely formatted fasta entry for the given sequence
sub fasta {
   my ($no) = @_;
   my $result = ">$sequenceheaders[$no]\n";
   for (my $i = 0; $i < length($sequences[$no]); $i += 100) {
      $result .= substr($sequences[$no], $i, 100);
      $result .= "\n"; }
   return $result;
}


# Sub for easy writing of html output
sub tag {
   my ($msg, @types) = @_;
   foreach my $type (@types) {
      if ($type eq 'p') {
         $msg .= "\n";
         $msg .= '<p>' if $www;
         $msg .= "\n"; }
      elsif ($type eq 'br') {
         $msg .= '<br>' if $www;
         $msg .= "\n"; }
      elsif ($type eq 'a') {
         $msg = "<a href=\"$msg\">$msg</a>" if $www; } 
      else {
         $msg = "<$type>$msg</$type>" if $www; } }
   return $msg;
}

sub version_copyright {
{print <<ENDPRINT
MaxAlign version 1.0
Increases the amount of information in an alignment by removing disturbing sequences

Copyright (C) 2007  Peter Wad Sackett, pws\@cbs.dtu.dk
                    Rodrigo de Oliveira, rodrigo\@cbs.dtu.dk

Official web site: http://www.cbs.dtu.dk/services/MaxAlign/

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
http://www.gnu.org/copyleft/gpl.html

ENDPRINT
}
exit;
}

### SUBS BELOW ARE NOT USED ###

# Calculates the alignment area for a given set the hard way
# Used for debugging purposes.
sub calcalignvol {
   my ($bits) =@_;
   my $freecols = 0;
   my $bst = unpack('b*', $bits);
   for (my $i = 0; $i < $longestseq; $i++) {
      my $gap = 0;
      for (my $j = 0; $j <= $#gapmatrix; $j++) {
         next if substr($bst, $j, 1) eq '1';
         next unless $gapmatrix[$j][$i];
         $gap = 1;
         last; }
      unless ($gap) {
         $freecols++;
         next; } }
   return $freecols * ($seqno - unpack('%32b*',$bits));
}

# Checks the hard way for a given set, that the gaps as bit vectors
# are correct. Used for debugging purposes.
sub check {
   my ($no) = @_;
   my $gaps = substr(unpack('b*', $gaps[$no]), 0, $longestseq);
   my $set = substr(unpack('b*', $sets[$no]), 0, $seqno);
   my $error = 0;
   for (my $g = 0; $g < $longestseq; $g++) {
#      die if $error;
      my $gapinset = 0;
      my $gapnotset = 0;
      for (my $i = 0; $i <= $seqno; $i++) {
         if ($gapmatrix[$i][$g]) {
	    if (substr($set, $i, 1) eq '1') {
	       $gapinset++; }
	    else {
	       $gapnotset++; } } }
      if (substr($gaps, $g, 1) eq '1') {
         if ($gapnotset) {
            print "gap is also OUTSIDE set $no, but should not be\n";
	    $error++; }
         unless ($gapinset) {
            print "gap is NOT INSIDE set $no, but should be\n";
	    $error++; } }
   }
}
      

# This simple heuristic greedily tries to calculate the best alignment area
# for a given problem without recalculation of the gapmatrix, sets and so
# forth, as the main heuristic does.
# It is not too good (and not really bad), but in any case it is not used.
sub simple_heuristic {
   my ($problem, $basesets, $basegaps) = @_;
   $problem = 'X' x scalar @sets unless defined $problem;
   $basegaps = pack('b*', '0' x $longestseq) unless defined $basegaps;
   $basesets = pack('b*', '0' x scalar @sets) unless defined $basesets;
   my $pointer = 0;
   while (-1 != ($pointer = index($problem, $pointer, '1'))) {
      $basegaps |= $gaps[$pointer];
      $basesets |= $sets[$pointer];
      $pointer++; }
   my $hmaxchar = ($freecolumns + unpack('%32b*', $basegaps)) * ($seqno - unpack('%32b*', $basesets));
   for (;;) {
      my %setorder;
      $pointer = 0;
      while (-1 != ($pointer = index($problem, 'X', $pointer))) {
         $setorder{$pointer} = ($freecolumns + unpack('%32b*', $basegaps | $gaps[$pointer])) *
                               ($seqno - unpack('%32b*', $basesets | $sets[$pointer]));
         $pointer++; }
      last unless %setorder;                               
      my $best = (sort {$setorder{$b} <=> $setorder{$a}} keys %setorder)[0];
      last if $hmaxchar >= $setorder{$best};
      $hmaxchar = $setorder{$best};
      $basegaps |= $gaps[$best];
      $basesets |= $sets[$best];
      substr($problem, $best, 1, '1'); }
   return $hmaxchar;
}


#### IDEAS TO BE TESTED ####

# If some (or one) sequences MUST be included in the optimizing, that will
# help a lot (depending on the sequences).

#### FAILED IDEAS ####
=pod
This is three different bound algorithms. None of them are really worthy of
the name, and they all failed.

      # First trial, just tests if there is something to improve. Don't use it
      my $canscoreincrease = 0;
      for (my $i = $pointer; $i < $setsno; $i++) {
         next unless substr($problem, $i, 1) eq 'X';
         if ($score < ($freecolumns + unpack('%32b*', $uniongaps | $gaps[$i])) *
                      ($seqno - unpack('%32b*', $unionsets | $sets[$i]))) {
            $canscoreincrease = 1;
	    last; } }
      if (not $canscoreincrease and $score < $alignmentarea) {
         $pruning++;
         if ($debug > 2) {
            print "Bound pruning; $problem\n"; }
         last; }

      # Second trial, add up positive improvements for individual sets, No good.
      my $addscore = 0;
      for (my $i = $pointer; $i < $setsno; $i++) {
         next unless substr($problem, $i, 1) eq 'X';
         my $newscore = ($freecolumns + unpack('%32b*', $uniongaps | $gaps[$i])) *
                        ($seqno - unpack('%32b*', $unionsets | $sets[$i]));
         if ($newscore > $score) {
            $addscore += $newscore - $score; } }
      if (not $addscore and ($score + $addscore) < $alignmentarea) {
         $pruning++;
         if ($debug>2) {
            print "Bound pruning; $problem\n"; }
         last; }

      # Third trial, pruning by number of iterations. Bad but quick
      if ($problem =~ tr/1/1/ > $maxremoveset and $score < $alignmentarea) {
         $pruning++;
         if ($debug>2 and $pruning % 10000 == 0) {
            print "Bound pruning; $problem\n"; }
         last; }
=cut


=pod
Here are a few failed ideas on how to eliminate sets.


# A solution is likely to remove around $iterationscount sets, since that is
# what the greedy heuristic did. Based on that we could eliminate "almost
# too large" sets, i.e. sets that in any combination with a number of other
# sets shows themselves to be too large.
# LATER COMMENT: This idea showed itself to be rubbish
$lastgapcolumns = -1;
while ($lastgapcolumns != $gapcolumns) {
   $lastgapcolumns = $gapcolumns;
   for (my $i = $#sets; $i >= 0; $i--) {
      my $almostlarge = 1;
      my $mainset = $sets[$i];
      for (my $j = $#sets; $j >= 0; $j--) {
         my $minorset = $sets[$j];
         my $union = unpack('b*', $mainset | $minorset);
         next if $union eq unpack('b*', $mainset) or $union eq unpack('b*', $minorset);
         next if $alignmentarea > ($longestseq - $gapcolumns) * ($seqno - $union =~ tr/1/1/);
         $almostlarge = 0;
         last; }
      if ($almostlarge) {
         splice(@gaps, $i, 1);
         splice(@sets, $i, 1);
         $gapcolumns = &getgapcolumns; } } }


# Are there any sets that MUST be removed in order to optimize the alignment
# score. These sets need not be considered, but simply eliminated and
# and new sets calculated.
# LATER COMMENT: This is a good idea and but it is hard figure out a safe
# way to do it. Are you REALLY REALLY sure that these sequences must be
# removed in a optimal solution? Perhaps the problem can be attacked from
# an other angle than this.
my $preunionsets = pack('b*', '0' x $seqno);
my $preuniongaps = pack('b*', '0' x $longestseq);
my $preproblem = 'X' x scalar @sets;
for (my $i = $#sets; $i >= 0; $i--) { #last;
   my $problem = 'X' x scalar @sets;
   substr($problem, $i, 1, '1');
   my $hremove = &heuristic($problem);
   substr($problem, $i, 1, '0');
   my $hstay = &heuristic($problem);
#   print "$hstay - $hremove\n";
   next if $hstay == $hremove;
   if ($hremove > $hstay) {
      next unless 0.9*$hremove > $hstay;
      substr($preproblem, $i, 1, '1');
      $preunionsets |= $sets[$i];
      $preuniongaps |= $gaps[$i]; 
      print "Union Set: ", unpack('%32b*', $preunionsets), "\n";
      print "Union Gap: ", unpack('%32b*', $preuniongaps), "\n";
      print "E MAXCHAR: ", ($seqno - unpack('%32b*', $preunionsets)) * ($freecolumns + unpack('%32b*', $preuniongaps)), "\n\n"; 
      }
   elsif (0) { #($hremove < $hstay) {
      splice(@gaps, $i, 1);
      splice(@sets, $i, 1);
      print "Eliminiation\n"; } 
   }
if ($debug>1) {
   print 'REM Heuristic: ', &heuristic($preproblem), "\n"; }
for (my $i = $#sets; $i >= 0; $i--) {
   if (unpack('b*', $preunionsets | $sets[$i]) eq unpack('b*', $preunionsets)) {
      $preuniongaps |= $gaps[$i];
      splice(@gaps, $i, 1);
      splice(@sets, $i, 1); } }
if ($debug>1) {
   print 'ELIM Heuristic: ', &heuristic(undef, $preunionsets, $preuniongaps), "\n"; 
   print "ELIM MAXCHAR: ", ($seqno - unpack('%32b*', $preunionsets)) * ($freecolumns + unpack('%32b*', $preuniongaps)), "\n"; 
   print "Sets after opportunistic elimination: ", scalar @sets, "\n"; }
&benchmark('Opportunistic Elimination');


# Some sets exhibits synergy, i.e. their removal is improving the alignment area
# more as a pair, than the improvement from their individual removal added together.
# LATER COMMENT: The idea is tempting, but wrong. Here it is proven for 2 sets.
# Assumption: The sets do NOT overlap. Calculation: What do we gain be removal.
# (SeqNo - SetSize1) * GapsInSet1 < (Length - GapsInSet1 - GapsInSet2) * SetSize1
# (SeqNo - SetSize2) * GapsInSet2 < (Length - GapsInSet1 - GapsInSet2) * SetSize2
# (SeqNo - SetSize1- SetSize2) * (GapsInSet1 + GapsInSet2) > (Length - GapsInSet1 - GapsInSet2) * (SetSize1 + SetSize2)
# Solving the inequalities gives: 0 <= -SetSize2 * GapsInSet1 - SetSize1 * GapsInSet2 : FALSE for any sets
# Assuming the sets DO overlap, then there is another more powerfull set that incorporates
# the best qualities of the 2 sets (unless one is a subset of another - then synergy is out)
# and will be chosen for removal before any of the 2 sets. When the powerfull set is removed
# then what (gaps) are left - none in common - i.e. we are back to the first assumption.
# Conclusion: Synergy is already built-in in the way the algorithm works - how the
# sets are calculated.
my %synergy;
my ($pos, $null, $neg) = (0, 0, 0);
%setorder = ();
for (my $i = $#sets; $i >= 0; $i--) {
   $setorder{$i} = ($seqno - unpack('%32b*', $sets[$i])) * ($freecolumns + unpack('%32b*', $gaps[$i])); }
for (my $i = $#sets; $i > 0; $i--) {
   for (my $j = $i-1; $j >= 0; $j--) {
      my $syn = ($seqno - unpack('%32b*', $sets[$i] & $sets[$j])) *
                ($freecolumns + unpack('%32b*', $gaps[$i] & $gaps[$j]));
      $syn = $setorder{$i}+$setorder{$j}-$orig_alignmentarea;
      $neg++ if $syn < 0;
      $null++ if $syn == 0;
      $pos++ if $syn > 0;
      $synergy{"$i, $j"} = $syn if $syn > 0; } }
print "Negative syn: $neg\n";
print "Null syn: $null\n";
print "Positive syn: $pos\n";
#foreach my $k (sort keys %synergy) {
#   print "$k: $synergy{$k}\n"; }
exit;
=cut