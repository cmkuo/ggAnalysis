# Script to do simple tests on source code level
# awk -f ntuplizer_checker.awk ggNtuplizer.cc

BEGIN{
	print "ggNtuple.cc sanity check";
	branchclear[0]=0;
	branchpush[0]=0;
}

/tree_->Branch/{
	#print "considering line "$0;
	if(match($0,/&.+\)/)) {
		bname = substr($0,RSTART+1,RLENGTH-2)
		print "adding branch "bname;
		branchclear[bname]=10;
		branchpush[bname]=10;
	}
}

/\.clear\(\)/{
	#print "considering line "$0;
	match($0,/[^ \t]+\./)
	bname = substr($0,RSTART,RLENGTH-1);
	#print bname
	if(branchclear[bname]) {
		print "cleared branch "bname;
		branchclear[bname]+=1;
	}
}

/\.push_back\(/{
	#print "considering line "$0;
	match($0,/\t*\ *\t*[^ \t]+\ *\t*\ *\.push_back\(/);
	bname1 = substr($0,RSTART,RLENGTH-11);
	#print bname1
	match(bname1,/[^ \t]+/);
	bname = substr(bname1,RSTART,RLENGTH);
	#print bname
	if(branchpush[bname]){
		print "push_back branch "bname;
		branchpush[bname]+=1;
	}
}
END{
	print " =================================   Report   ======================================= "
	for(bname in branchpush){
                if(branchpush[bname] && branchpush[bname] > 11) print "Branch "bname" was push_back-ed more that once : "branchpush[bname]-10
        }
	print "======================================================================================"
	for(bname in branchclear){
		if(branchclear[bname] && branchpush[bname]>10 && branchclear[bname] == 10) print bname".clear();"
	}
	print "need these lines ^^^=================================================================="

	for(bname in branchclear){
		if(branchclear[bname] && branchpush[bname]>10 && branchclear[bname] > 11) print bname" was cleared more than once!!!"
	}

}

