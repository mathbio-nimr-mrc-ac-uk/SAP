#!/usr/bin/perl

#==============================================================================
# $Id$
# sap.example.pl : run SAP WEB server application
# (C) 2008 Jens Kleinjung and Alessandro Pandini
#==============================================================================

$phpdir = "/srv/www/htdocs/wiki/php";

#______________________________________________________________________________
# copy default files to working dir

system("cp /tmp/databases/examples/sap.parameters parameters");
$pwddir = $ENV{'PWD'};                                                                                                            
@pwdpath = split('/',$pwddir);
$jobId = $pwdpath[-1];
system("echo mbjob.uniqid=$jobId >> parameters");

#______________________________________________________________________________
# read WEB parameters
open(IN, "parameters");
@pars = <IN>;
close(IN);

foreach $par (@pars)
{
    chomp $par;
    ($tag, $val) = split(/\=/, $par);
    # pdbin1
    if ($tag eq "file.PDBINPUT1.name")
    { $pdbin1 = $val;
      $basename1 = substr($pdbin1, 0, $#pdbin1 - 3);
    };
    # pdbin2
    if ($tag eq "file.PDBINPUT2.name")
    { $pdbin2 = $val;
      $basename2 = substr($pdbin2, 0, $#pdbin2 - 3);
    };
    # PDBID1
    if ($tag eq "proteinID1")
    { $proteinID1 = $val};
    # PDBID2
    if ($tag eq "proteinID2")
    { $proteinID2 = $val};
    # job ID
    if ($tag eq "mbjob.uniqid")
    { $jobId = $val};
    # user email address
    if ($tag eq "useremail")
        { $useremail = $val};
}

#______________________________________________________________________________
# run pdbEncode application

# compose command string
$command_string = "$phpdir/binary/sap";

if ($pdbin1 ne ""){
    $command_string .= " $pdbin1";
}else{
	&process_pdb($proteinID1);
    $command_string .= " $proteinID1.pdb";
}

if ($pdbin2 ne ""){
    $command_string .= " $pdbin2";
}else{
	&process_pdb($proteinID2);
    $command_string .= " $proteinID2.pdb";
}

$command_string .= " > aln.log";

# execute command string
system($command_string);

# flag up status
system("touch job.completed");
system("echo '0' > status");

&send_email();

# subroutine: send email
sub send_email
{
        system("echo 'SAP: Job ID $jobId completed. Go to http://mathbio.nimr.mrc.ac.uk/wiki/SAP, insert the job ID and press -Show Results- to retrieve your results.' | mail -s 'SAP job completion' $useremail");
        system("echo 'SAP: Job ID $jobId completed. Go to http://mathbio.nimr.mrc.ac.uk/wiki/SAP, insert the job ID and press -Show Results- to retrieve your results.' > mail_sent.txt");
}

# load PDB file from PDB resource
sub process_pdb
{
	$query = $_[0];
    # link to pdb structure and unzip
    $pdbquery = "$query" . ".pdb.gz";
    $unzipquery = "$query" . ".pdb";
    print("Retrieving $pdbquery\n");
    $fileUrl = "http://www.rcsb.org/pdb/files/" . $pdbquery;
    $saveAs = "./" . $pdbquery;
    print("wget -v -O $saveAs $fileUrl\n");
    $error = system("export http_proxy=http://proxy.nimr.mrc.ac.uk:8000; wget -v -O $saveAs $fileUrl");
	print("Download error status: $error\n");
    if ($error > 0) {
		print("Warning: Failed to download $query!\n");
		system("echo '2' > status");
		exit(1);
    }
    system("gunzip $pdbquery");
}
