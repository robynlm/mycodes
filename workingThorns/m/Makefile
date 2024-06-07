# -*-Makefile-*-

Analytic.out: Analytic.m
	\rm -rf CT_Analytic
	./runmath.sh $^
	./tl.patch
	for thorn in CT_Analytic; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done

Dust.out: Dust.m
	\rm -rf CT_Dust
	./runmath.sh $^
	./dust.patch
	for thorn in CT_Dust; do				\
		./copy-if-changed.sh $$thorn ../$$thorn;	\
	done
	echo "copy-if-changed was run"
