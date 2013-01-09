CXX=g++
%.o: seven-op.h analysis.h
seven-op: seven-op.o seven-op-common.o analysis.o
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LDFLAGS)
#coordinate only
%.ar3a: %.gro
	./seven-op -c < $< > $@
#network topology only
%.ngph: %.gro
	./seven-op -n < $< > $@
#order parameters, cluster info, and their spatial correlations
%.op: %.gro
	./seven-op -r < $< > $@
#raw order parameter only
%.op0: %.gro
	./seven-op -0 < $< > $@
#network difference from the given file (0ps.ngph)
%.diff: %.ngph %.ar3a 0ps.ngph
	cat $*.ngph $*.ar3a | python netdiff.py 0ps.ngph > $@
#spatial correlation of the order parameter
%.op.scf: %.op
	awk '(NF==4){x=int($$3*3+0.5);s[x]+=$$4;n[x]++}END{for(i in s){print i/3.0, s[i]/n[i], n[i]}}' $< | sort -n > $@
%.op0.tcf: %.op0 0ps.op0
	paste $^ | tail +3 | awk '{s+=$$1*$$2;n++}END{print s/n}' > $@


