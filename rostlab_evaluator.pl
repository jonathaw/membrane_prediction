#!/usr/bin/perl
use strict;
use List::Util qw[max];

#################
###Subroutines###
#################
sub getPerSegmentBaseMeasures;
sub getHelixBegin;
sub getHelixEnd;
sub getHelixLength;
sub getNumHelices;

################
###Variables####
################
##Constants
my $MARK_S_HELIX = 'H|h|C|L';
my %reverseHelixSymbols = (
    'H' => 'h|C|L', 'h' => 'H|C|L',
    'C' => 'h|H|L', 'L' => 'H|h|C'
);


#Being/end residues of predicted may differ by at most $MAX_END_OFFSET from the observed helix
my $MAX_END_OFFSET = 5;

#Read obs and pred from command line (pipe)
my $obs =$ARGV[0];
my $pred =$ARGV[1];
#Example 1
#my $obs ='1111HHHHHHHHHHHHHHH111111111111HHHHHHHHHHHHHH11111111111111111111111111111111111111111111HHHHHHHHHHHHH111111111111111HHHHHHHHHHHHHH111111111HHHHHHHHHHHHH111111111111111111HHHHHHHHHHHHH111111111111111111111111111111111111111111111HHHHHHHHHHHHHHHHHHH111111HHHHHHHHHHHHHHH11111111111111111111111HHHHHHHHHHHHHH11111111HHHHHHHHHHHHHH1111111111111HHHHHHHHHHHHHHH1111111111111111111111111111111111HHHHHHHHHHHHHHHHHH11111111111111HHHHHHHHHHHHHH1111111111111111111111111111111HHHHHHHHHHHHHHH111111111111111111111111111';
#my $pred = "1HHHHHHHHHHHHHHHHHHH111111111HHHHHHHHHHHHHHHHHHHH111111111111111111111111111111111111HHHHHHHHHHHHHHHHHHHHH11111111HHHHHHHHHHHHHHHHHHH1111HHHHHHHHHHHHHHHHHHHH1111111111111111HHHHHHHHHHHHHHHHHHHHHHHHH11111111111111111111HHHHHHHHHHHHHHHHHHHHHH111111111111HHHHHHHHHHHHHHHHHHHH111111111111HHHHHHHHHHHHHHHHHHHH1111111111HHHHHHHHHHHHHHHHHHHH1111111HHHHHHHHHHHHHHHHHHHHHH1111111111111111111HHHHHHHHHHHHHHHHHHHHHHH11111111111HHHHHHHHHHHHHHHHHHHHHHH1111111111111111111111111HHHHHHHHHHHHHHHHHHHH1111111111111111111111111";

#Example 2: lower case 'h' is used, when an annotated helix directly follows another to distinguish them as two separate helices
# my $obs  = "1111HHHHHHHHHhhhhhhhhhhhhhhHHHHHHHH111111";
# my $pred = "HHHHHHHHHHH1HHHHHHHHHHHHH1111HHHHHHHHHH11";

my ($corrPred) = getPerSegmentBaseMeasures($obs, $pred);

print "Number of correctly predicted helices: $corrPred\n\n";

###################
###Subroutines#####
###################
sub getPerSegmentBaseMeasures
{
    my $obs = shift;
    my $pred = shift;


    my @split_obs = split(//, $obs);
    my @split_pred = split(//, $pred);
    my ($obs_c, $pred_c);
    my $numCorrPred = 0;

    for (my $i = 0; $i < scalar(@split_obs); $i++)
    {
        $obs_c = $split_obs[$i];
        $pred_c = $split_pred[$i];

        if($obs_c =~ m/$MARK_S_HELIX/)
        {
            if($pred_c =~ m/$MARK_S_HELIX/)
            {
                if(&checkSufficientOverlap($obs, $pred, $i))
                {
                    $numCorrPred++;

                    #We have found two helices with sufficient overlap, the observed helix has been counted as correctly predicted
                    #Now advance the index to the next relevant point
                    #Since the current observed helix has been handled we can move at least to the end index of this helix
                    #Unless the current predicted helix is longer, then move to the end index of that one, since one predicted helix, can never 'count' for several observed ones (the end point constraint would take care of this anyway)
                    $i = max(getHelixEnd($obs, $i), getHelixEnd($pred, $i));    #No +1 here, since, the loop is not over and will auto increment $i anyway
                }
            }
        }
    }
    return ($numCorrPred);
}

sub checkSufficientOverlap
{
    my $obs = shift;
    my $pred = shift;
    my $ind = shift;
#     my $currentOverlap = shift;

    #By default, assume sufficiency
    my $suff = 1;

    #Get begin and end indices of helices at current index
    my $obsbeg = getHelixBegin($obs, $ind);
    my $obsend = getHelixEnd($obs, $ind);

    my $predbeg = getHelixBegin($pred, $ind);
    my $predend = getHelixEnd($pred, $ind);

    #Condition 1: Sufficient overlap, not used

    #Condition 2: Beginning and end of both helices must not deviate by more than $MAX_END_OFFSET
    #This is implemented so that $MAX_END_OFFSET is the number of residues from the one first/last residue, to the other first/last residue, both inclusive
    #Example: For $MAX_END_OFFSET = 3, this is the maximum that is still ok, because abs(0 - 2) +1 = 3:
    #PRD:CCH
    #OBS:HHH
    if(abs($obsbeg-$predbeg)+1 > $MAX_END_OFFSET || abs($obsend - $predend)+1 > $MAX_END_OFFSET)
    {
        #Too large deviation in end points
        $suff = 0;
    }

    return $suff;
}

sub getHelixEnd
{
    my $struct = shift;
    my $index = shift; #Some index of the helix for which we need the end

    my $end = $index;
    my $helixMarker = '';
    my @split_struct = split(//, $struct);
    for(my $i=$index; $i < scalar(@split_struct); $i++)
    {
        if($helixMarker eq '' && $split_struct[$i] =~ m/$MARK_S_HELIX/)
        {
            #The helix we want to look at is here, so save the marker (either H or h) that is used for this specific helix
            $end = $i;
            $helixMarker = $split_struct[$i];
        }
        else
        {
            if($helixMarker ne '' && $split_struct[$i] eq $helixMarker)
            {
                #Since we are not at the first index we look at anymore, we only allow the marker that is used for this specific helix from now on
                $end = $i;
            }
            else
            {
                last;
            }
        }
    }

    return $end;
}

sub getHelixBegin
{
    my $struct = shift;
    my $index = shift; #Some index of the helix for which we need the beginning

    my $begin = $index;
    my @split_struct = split(//, $struct);
    my $helixMarker = '';
    for(my $i=$index; $i >= 0; $i--)
    {
        if($helixMarker eq '' && $split_struct[$i] =~ m/$MARK_S_HELIX/)
        {
            #The helix we want to look at is here, so save the marker (either H or h) that is used for this specific helix
            $begin = $i;
            $helixMarker = $split_struct[$i];
        }
        else
        {
            if($helixMarker ne '' && $split_struct[$i] eq $helixMarker)
            {
                #Since we are not at the first index we look at anymore, we only allow the marker that is used for this specific helix from now on
                $begin = $i;
            }
            else
            {
                last;
            }
        }
    }

    return $begin;
}

sub getNumHelices
{
    my $struct = shift;

    my @split_struct = split(//, $struct);
    my $hn = 0;
    my $intmh = 0;
    my $currHelixSymb = '';

    for ( my $loop = 0; $loop < scalar( (@split_struct) ); $loop++ )
    {
        my $ss = $split_struct[$loop];

        if ( $ss =~ m/$MARK_S_HELIX/ )
        {
            if ( !$intmh )
            {
                $hn++;
                $currHelixSymb = $ss;
            }
            else
            {
                my $othersymb = $reverseHelixSymbols{$currHelixSymb};
                if($ss =~ m/$othersymb/)
                {
                    $currHelixSymb = $ss;
                    $hn++;
                }
            }

            $intmh = 1;
        }
        else
        {
            $intmh = 0;
            $currHelixSymb = '';
        }
    }

    return $hn;
}


1;