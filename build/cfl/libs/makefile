#     $Id$
#=============================================================================
#
#              creates common libraries for subroutines common to 
#              both the sequential and parallel versions of cfl3d
#              note that this excludes any file that contains 
#              MPI/sequential protocols activated by "if def"
#
#=============================================================================

# ****************************** CREATE LINKS ********************************

link: lncode lnhead

lncode:
	@ echo "        linking source code"
	ln -sf $(CFLSRC_S)/*.F .
	ln -sf $(CFLSRC_S)/*.F90 .
	ln -sf $(CFLSRC_S)/*.c .

lnhead:

# ****************************** SUFFIX RULES ********************************
.SUFFIXES:
.SUFFIXES:.F .o .F90 .c
.F.o: 
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F
.F90.o: 
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F90
.c.o:
	$(CC) $(CFLAG) -c $*.c

# *************************** LIBRARY DEFINITIONS ****************************

FSRC_LIBS = \
abciz.F                \
abcjz.F                \
abckz.F                \
add2x.F                \
aesurf.F               \
ae_pred.F              \
af3f.F                 \
amafi.F                \
amafj.F                \
amafk.F                \
arc.F                  \
arclen.F               \
augmntq.F              \
avghole.F              \
avgint.F               \
barth3d.F              \
bc.F                   \
bc1000.F               \
bc1001.F               \
bc1002.F               \
bc1003.F               \
bc1005.F               \
bc1008.F               \
bc1011.F               \
bc1012.F               \
bc1013.F               \
bc2002.F               \
bc2003.F               \
bc2004.F               \
bc2005.F               \
bc2005i_d.F            \
bc2005j_d.F            \
bc2005k_d.F            \
bc2006.F               \
bc2007.F               \
bc2008.F               \
bc2009.F               \
bc2010.F               \
bc2016.F               \
bc2019.F               \
bc2026.F               \
bc2102.F               \
bc9999.F               \
bcchk.F                \
bcnonin.F              \
bc_delt.F              \
bc_info.F              \
bc_vdsp.F              \
bc_xmera.F             \
blkmax.F               \
blnkfr.F               \
blocki.F               \
blocki_d.F             \
blockj.F               \
blockj_d.F             \
blockk.F               \
blockk_d.F             \
blomax.F               \
bsub.F                 \
bsubp.F                \
cblki.F                \
cblki_d.F              \
cblkj.F                \
cblkj_d.F              \
cblkk.F                \
cblkk_d.F              \
ccf.F                  \
ccomplex.F             \
cctogp.F               \
cellvol.F              \
cgnstools.F            \
chkdef.F               \
chkrap.F               \
chkrot.F               \
chkroti_d.F            \
chkrotj_d.F            \
chkrotk_d.F            \
chksym.F               \
cntsurf.F              \
coll2q.F               \
collapse.F             \
colldat.F              \
collmod.F              \
collq.F                \
collqc0.F              \
collv.F                \
collx.F                \
collxt.F               \
collxtb.F              \
csout.F                \
csurf.F                \
ctime1.F               \
dabciz.F               \
dabcjz.F               \
dabckz.F               \
deform.F               \
delintr.F              \
delq.F                 \
delv.F                 \
delv_k.F               \
delv_p.F               \
dfbtr.F                \
dfbtrp.F               \
dfhat.F                \
dfluxpm.F              \
diagi.F                \
diagj.F                \
diagk.F                \
diagnos.F              \
dird.F                 \
direct.F               \
dlutr.F                \
dlutrp.F               \
dsmin.F                \
dthole.F               \
expand.F               \
extra.F                \
extrae.F               \
fa.F                   \
fa2xi.F                \
fa2xj.F                \
fa2xk.F                \
fcd.F                  \
fdelay.F               \
ffluxl.F               \
ffluxr.F               \
ffluxr_cd.F            \
ffluxv.F               \
ffluxv1.F              \
fhat.F                 \
fill.F                 \
fluxm.F                \
fluxp.F                \
fmaps.F                \
force.F                \
foureqn.F              \
genforce.F             \
getdelt.F              \
getdhdr.F              \
getibk.F               \
getibk0.F              \
getsurf.F              \
get_bvals.F            \
gfluxl.F               \
gfluxr.F               \
gfluxr_cd.F            \
gfluxv.F               \
gfluxv1.F              \
global.F               \
global0.F              \
global2.F              \
gradinfo.F             \
grdmove.F              \
hfluxl.F               \
hfluxr.F               \
hfluxr_cd.F            \
hfluxv.F               \
hfluxv1.F              \
histout.F              \
hole.F                 \
i2x.F                  \
i2xi_d.F               \
i2xj_d.F               \
i2xk_d.F               \
i2xs.F                 \
i2xsi_d.F              \
i2xsj_d.F              \
i2xsk_d.F              \
init.F                 \
initnonin.F            \
initvist.F             \
init_ae.F              \
init_mast.F            \
init_rb.F              \
init_trim.F            \
int2.F                 \
int2_d.F               \
intrbc.F               \
invert.F               \
l2norm.F               \
l2norm2.F              \
lamfix.F               \
ld_dati.F              \
ld_datj.F              \
ld_datk.F              \
ld_qc.F                \
lead.F                 \
lesdiag.F              \
loadgr.F               \
metric.F               \
mms.F                  \
moddefl.F              \
modread.F              \
mreal.F                \
mvdat.F                \
my_flush.F             \
newfit.F               \
outbuf.F               \
parser.F               \
pltmode.F              \
pre_blockbc.F          \
pre_blocki.F           \
pre_blockj.F           \
pre_blockk.F           \
pre_cblki.F            \
pre_cblkj.F            \
pre_cblkk.F            \
pre_embed.F            \
pre_patch.F            \
pre_period.F           \
project.F              \
prolim.F               \
prolim2.F              \
q8sdot.F               \
q8smax.F               \
q8smin.F               \
q8vrev.F               \
qface.F                \
rb_corr.F              \
rb_pred.F              \
rcfl.F                 \
readdat.F              \
readkey.F              \
rechk.F                \
resadd.F               \
resid.F                \
resnonin.F             \
rie1d.F                \
rie1de.F               \
rotate.F               \
rotateq.F              \
rotateq0.F             \
rotateq2_d.F           \
rotateqb.F             \
rotatmc.F              \
rotatp.F               \
rotsurf.F              \
rp3d.F                 \
rpatch.F               \
rpatch0.F              \
rsmooth.F              \
rsurf.F                \
setblk.F               \
setcorner.F            \
setdqc0.F              \
setqc0.F               \
setseg.F               \
shear.F                \
sijrate2d.F            \
sijrate3d.F            \
spalart.F              \
swafi.F                \
swafj.F                \
swafk.F                \
tau.F                  \
tau2x.F                \
tdq.F                  \
tfiedge.F              \
tfiface.F              \
tfivol.F               \
threeeqn.F             \
tinvr.F                \
tmetric.F              \
topol.F                \
topol2.F               \
trace.F                \
trans.F                \
transmc.F              \
transp.F               \
triv.F                 \
trnsurf.F              \
twoeqn.F               \
unld_qc.F              \
update.F               \
u_doubleprime.F        \
varnam.F               \
vargrad.F              \
vlutr.F                \
vlutrp.F               \
wkstn.F                \
wmag.F                 \
xe.F                   \
xe2.F                  \
xlim.F                 \
xlsfree.F              \
xmukin.F               \
xtbatb.F               \
xupdt.F                \
xyzintr.F              

F90SRC_LIBS = module_profileout.F90 module_contour.F90  module_kwstm.F90 module_stm_2005.F90

FSRC_SPEC = addx.F

CSRC_LIBS = bessel.c

FOBJ_LIBS = $(FSRC_LIBS:.F=.o)

FOBJ_SPEC = $(FSRC_SPEC:.F=.o)

COBJ_LIBS = $(CSRC_LIBS:.c=.o)

F90OBJ_LIBS = $(F90SRC_LIBS:.F90=.o) 

COMMONLIB = libcommon.a

initvist.o: module_kwstm.o
resid.o: module_kwstm.o
init.o:  module_kwstm.o
module_kwstm.o: module_profileout.o
bc2005.o: module_stm_2005.o
bc2005i_d.o: module_stm_2005.o
bc2005j_d.o: module_stm_2005.o
bc2005k_d.o: module_stm_2005.o

$(COMMONLIB): $(FSRC_LIBS) $(FOBJ_LIBS) $(FSRC_SPEC) $(FOBJ_SPEC) \
	$(CSRC_LIBS) $(COBJ_LIBS) $(F90OBJ_LIBS)
	ar $(AROPT) $(COMMONLIB) $(FOBJ_LIBS) $(F90OBJ_LIBS) $(FOBJ_SPEC) $(COBJ_LIBS) 
	@$(RANLIB) $(COMMONLIB)

$(FOBJ_LIBS): 
	$(FTN) $(CPPOPT) $(FFLAG) -c $*.F

$(FOBJ_SPEC):
	$(FTN) $(CPPOPT) $(FFLAG_SPEC) -c $*.F

$(COBJ_LIBS):
	$(CC) $(CFLAG) -c $*.c

$(EXEC): $(COMMONLIB)
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "
	@ echo "               DONE:  $(COMMONLIB) created                    "
	@ echo "                                                              "
	@ echo "           the archive libraries can be found in:             "
	@ echo "                                                              "
	@ echo "                   $(DIR)/$(COMMONLIB)                        "
	@ echo "                                                              "
	@ echo "=============================================================="
	@ echo "                                                              "

# ****************************** CLEAN/SCRUB *********************************

# the @touch is used to (silently) create some temp files to prevent irksome
# warning messages are sometimes created if there are no *.whatever files and
# one tries to remove them

cleano:
	@touch temp.o
	-rm -f *.o *.mod

cleane:
	-rm -f $(EXEC) 

cleana:
	@touch temp.a
	-rm -f *.a

cleanf:
	@touch temp.f
	-rm -f *.f

cleanh:
	@touch temp.h
	-rm -f *.h

cleang:
	@touch temp.F
	-rm -f *.F *.F90 

cleanc:
	@touch temp.c
	-rm -f *.c

scrub: cleana cleano cleane cleanf cleanh cleang cleanc
