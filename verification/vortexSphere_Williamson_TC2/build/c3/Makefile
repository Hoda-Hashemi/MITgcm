# Multithreaded + multi-processing makefile for:
#    Linux anode07 3.10.0-1062.18.1.el7.x86_64 #1 SMP Tue Mar 17 23:49:17 UTC 2020 x86_64 x86_64 x86_64 GNU/Linux
# This makefile was generated automatically on
#    Thu Jun  4 17:46:47 EEST 2026
# by the command:
#    /home/hmh85/scratch/MITgcm/tools/genmake2 -rootdir=/home/hmh85/scratch/MITgcm -mpi -mods=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code -of=/home/hmh85/scratch/MITgcm/tools/build_options/linux_amd64_gfortran
# executed by:
#    hmh85@anode07:/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/build/c3
#
# OPTFILE="/home/hmh85/scratch/MITgcm/tools/build_options/linux_amd64_gfortran"
#
# BUILDDIR     : Directory where object files are written
# SOURCEDIRS   : Directories containing the source (.F) files
# INCLUDEDIRS  : Directories containing the header-source (.h) files
# EXEDIR       : Directory where executable that is generated is written
# EXECUTABLE   : Full path of executable binary
#
# CPP          : C-preprocessor command
# INCLUDES     : Directories searched for header files
# DEFINES      : Macro definitions for CPP
# MAKEDEPEND   : Dependency generator
# FC           : Fortran compiler command
# FFLAGS       : Configuration/debugging options for FC
# FOPTIM       : Optimization options for FC
# LINK         : Command for link editor program
# LIBS         : Library flags /or/ additional optimization/debugging flags

ROOTDIR     = /home/hmh85/scratch/MITgcm
BUILDDIR    = .
SOURCEDIRS  =  /home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code $(ROOTDIR)/pkg/debug $(ROOTDIR)/pkg/diagnostics $(ROOTDIR)/pkg/generic_advdiff $(ROOTDIR)/pkg/mdsio $(ROOTDIR)/pkg/mom_common $(ROOTDIR)/pkg/mom_fluxform $(ROOTDIR)/pkg/mom_vecinv $(ROOTDIR)/pkg/monitor $(ROOTDIR)/pkg/mypackage $(ROOTDIR)/pkg/rw $(ROOTDIR)/eesupp/src $(ROOTDIR)/model/src
INCLUDEDIRS =  /apps/sw/mpi/openmpi/4.1.4-gcc-9.1.0-slurm-18.08.6/include /home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code $(ROOTDIR)/pkg/debug $(ROOTDIR)/pkg/diagnostics $(ROOTDIR)/pkg/generic_advdiff $(ROOTDIR)/pkg/mdsio $(ROOTDIR)/pkg/mom_common $(ROOTDIR)/pkg/mom_fluxform $(ROOTDIR)/pkg/mom_vecinv $(ROOTDIR)/pkg/monitor $(ROOTDIR)/pkg/mypackage $(ROOTDIR)/pkg/rw $(ROOTDIR)/eesupp/inc $(ROOTDIR)/model/inc
EXEDIR      = .
EXECUTABLE  = $(EXEDIR)/mitgcmuv
TOOLSDIR    = $(ROOTDIR)/tools
OADTOOLS    = 

#eh3  new defines for the adjoint work
AUTODIFF    = $(ROOTDIR)/pkg/autodiff
EXE_ADJ     = mitgcmuv_ad
EXE_TLM     = mitgcmuv_ftl
EXE_SVD     = mitgcmuv_svd

ENABLED_PACKAGES =  -DALLOW_DEBUG -DALLOW_DIAGNOSTICS -DALLOW_GENERIC_ADVDIFF -DALLOW_MDSIO -DALLOW_MOM_COMMON -DALLOW_MOM_FLUXFORM -DALLOW_MOM_VECINV -DALLOW_MONITOR -DALLOW_MYPACKAGE -DALLOW_RW
DISABLED_PACKAGES =  -UALLOW_ADMTLM -UALLOW_AIM_V23 -UALLOW_ATM2D -UALLOW_ATM_COMMON -UALLOW_ATM_COMPON_INTERF -UALLOW_ATM_OCN_COUPLER -UALLOW_ATM_PHYS -UALLOW_AUTODIFF -UALLOW_BBL -UALLOW_BLING -UALLOW_BULK_FORCE -UALLOW_CAL -UALLOW_CD_CODE -UALLOW_CFC -UALLOW_CHEAPAML -UALLOW_CHRONOS -UALLOW_COMPON_COMMUNIC -UALLOW_COST -UALLOW_CTRL -UALLOW_DIC -UALLOW_DOWN_SLOPE -UALLOW_EBM -UALLOW_ECCO -UALLOW_EMBED_FILES -UALLOW_EXCH2 -UALLOW_EXF -UALLOW_FIZHI -UALLOW_FLT -UALLOW_FRAZIL -UALLOW_GCHEM -UALLOW_GGL90 -UALLOW_GMREDI -UALLOW_GRDCHK -UALLOW_GRIDALT -UALLOW_ICEFRONT -UALLOW_KL10 -UALLOW_KPP -UALLOW_LAND -UALLOW_LAYERS -UALLOW_LONGSTEP -UALLOW_MATRIX -UALLOW_MNC -UALLOW_MY82 -UALLOW_OBCS -UALLOW_OBSFIT -UALLOW_OCN_COMPON_INTERF -UALLOW_OFFLINE -UALLOW_OPENAD -UALLOW_OPPS -UALLOW_PP81 -UALLOW_PROFILES -UALLOW_PTRACERS -UALLOW_RBCS -UALLOW_REGRID -UALLOW_RUNCLOCK -UALLOW_SALT_PLUME -UALLOW_SBO -UALLOW_SEAICE -UALLOW_SHAP_FILT -UALLOW_SHELFICE -UALLOW_SHOWFLOPS -UALLOW_SMOOTH -UALLOW_SPHERE -UALLOW_STEEP_ICECAVITY -UALLOW_STREAMICE -UALLOW_TAPENADE -UALLOW_THSICE -UALLOW_ZONAL_FILT

# These files are created by Makefile
SPECIAL_FILES = PACKAGES_CONFIG.h AD_CONFIG.h FC_NAMEMANGLE.h BUILD_INFO.h
EMBEDDED_FILES = 
# Unix ln (link)
LN = ln -s
# Dependency generator
MAKEDEPEND = $(TOOLSDIR)/xmakedepend
# Special preprocessor (KAP on DECs, FPP on Crays)
KPP = 
# Fortran compiler
FC = mpif77
# Fortran compiler
F90C = mpif90
# C compiler
CC = mpicc
# Link editor
LINK = mpif77 

# Defines for CPP
DEFINES = -DWORDLENGTH=4 -DNML_TERMINATOR -DALLOW_USE_MPI -DHAVE_SYSTEM -DHAVE_FDATE -DHAVE_ETIME_SBR -DHAVE_CLOC -DHAVE_SETRLSTK -DHAVE_SIGREG -DHAVE_STAT -DHAVE_FLUSH 
# Includes for CPP
INCLUDES =  -I/apps/sw/mpi/openmpi/4.1.4-gcc-9.1.0-slurm-18.08.6/include
# default CPP includes for makedepend
CPPINCLUDES =  -I/apps/sw/gcc-10.1.0/lib/gcc/x86_64-pc-linux-gnu/10.1.0/include -I/usr/local/include -I/apps/sw/gcc-10.1.0/include -I/apps/sw/gcc-10.1.0/lib/gcc/x86_64-pc-linux-gnu/10.1.0/include-fixed -I/usr/include
# Flags for KPP
KFLAGS1 = 
KFLAGS2 = 
# Optim./debug for FC
FFLAGS =  -fallow-argument-mismatch -fconvert=big-endian -fimplicit-none -mcmodel=medium
FOPTIM = -O3 -funroll-loops 
# Optim./debug for FC
F90FLAGS =  -fallow-argument-mismatch -fconvert=big-endian -fimplicit-none -mcmodel=medium
F90OPTIM = -O3 -funroll-loops
F90FIXEDFORMAT = -ffixed-form
# Flags for CC
CFLAGS = -O0 -mcmodel=medium
# Files that should not be optimized
NOOPTFILES =  ini_masks_etc.F
NOOPTFLAGS = -O0
# Flags and libraries needed for linking
LIBS = 
# Name of the makefile
MAKEFILE = Makefile
# Name of the make program
MAKE = make

F77_SRC_FILES =  \
 chksum_tiled.F \
 debug_call.F \
 debug_cs_corner_uv.F \
 debug_enter.F \
 debug_fld_stats_rl.F \
 debug_fld_stats_rs.F \
 debug_leave.F \
 debug_msg.F \
 debug_stats_rl.F \
 debug_stats_rs.F \
 fill_in_corners_rl.F \
 write_fullarray_rl.F \
 write_fullarray_rs.F \
 diag_calc_psivel.F \
 diag_cg2d.F \
 diag_vegtile_fill.F \
 diagnostics_addtolist.F \
 diagnostics_calc_phivel.F \
 diagnostics_check.F \
 diagnostics_clear.F \
 diagnostics_fill.F \
 diagnostics_fill_field.F \
 diagnostics_fill_rs.F \
 diagnostics_fill_state.F \
 diagnostics_fract_fill.F \
 diagnostics_ini_io.F \
 diagnostics_init_early.F \
 diagnostics_init_fixed.F \
 diagnostics_init_varia.F \
 diagnostics_interp_p2p.F \
 diagnostics_interp_vert.F \
 diagnostics_is_on.F \
 diagnostics_list_check.F \
 diagnostics_main_init.F \
 diagnostics_mnc_out.F \
 diagnostics_out.F \
 diagnostics_read_pickup.F \
 diagnostics_readparms.F \
 diagnostics_scale_fill.F \
 diagnostics_scale_fill_rs.F \
 diagnostics_set_calc.F \
 diagnostics_set_levels.F \
 diagnostics_set_pointers.F \
 diagnostics_setdiag.F \
 diagnostics_status_error.F \
 diagnostics_sum_levels.F \
 diagnostics_summary.F \
 diagnostics_switch_onoff.F \
 diagnostics_utils.F \
 diagnostics_write.F \
 diagnostics_write_adj.F \
 diagnostics_write_pickup.F \
 diagstats_ascii_out.F \
 diagstats_calc.F \
 diagstats_clear.F \
 diagstats_close_io.F \
 diagstats_fill.F \
 diagstats_global.F \
 diagstats_ini_io.F \
 diagstats_local.F \
 diagstats_mnc_out.F \
 diagstats_others_calc.F \
 diagstats_output.F \
 diagstats_set_pointers.F \
 diagstats_set_regions.F \
 diagstats_setdiag.F \
 gad_advection.F \
 gad_advscheme.F \
 gad_biharm_r.F \
 gad_biharm_x.F \
 gad_biharm_y.F \
 gad_c2_adv_r.F \
 gad_c2_adv_x.F \
 gad_c2_adv_y.F \
 gad_c2_impl_r.F \
 gad_c4_adv_r.F \
 gad_c4_adv_x.F \
 gad_c4_adv_y.F \
 gad_calc_rhs.F \
 gad_check.F \
 gad_del2.F \
 gad_diagnostics_init.F \
 gad_diagnostics_state.F \
 gad_diff_r.F \
 gad_diff_x.F \
 gad_diff_y.F \
 gad_dst2u1_adv_r.F \
 gad_dst2u1_adv_x.F \
 gad_dst2u1_adv_y.F \
 gad_dst2u1_impl_r.F \
 gad_dst3_adv_r.F \
 gad_dst3_adv_x.F \
 gad_dst3_adv_y.F \
 gad_dst3fl_adv_r.F \
 gad_dst3fl_adv_x.F \
 gad_dst3fl_adv_y.F \
 gad_dst3fl_impl_r.F \
 gad_exch_som.F \
 gad_fluxlimit_adv_r.F \
 gad_fluxlimit_adv_x.F \
 gad_fluxlimit_adv_y.F \
 gad_fluxlimit_impl_r.F \
 gad_grad_x.F \
 gad_grad_y.F \
 gad_implicit_r.F \
 gad_init_fixed.F \
 gad_init_varia.F \
 gad_os7mp_adv_r.F \
 gad_os7mp_adv_x.F \
 gad_os7mp_adv_y.F \
 gad_osc_hat_r.F \
 gad_osc_hat_x.F \
 gad_osc_hat_y.F \
 gad_osc_mul_r.F \
 gad_osc_mul_x.F \
 gad_osc_mul_y.F \
 gad_plm_fun.F \
 gad_ppm_adv_r.F \
 gad_ppm_adv_x.F \
 gad_ppm_adv_y.F \
 gad_ppm_flx_r.F \
 gad_ppm_flx_x.F \
 gad_ppm_flx_y.F \
 gad_ppm_fun.F \
 gad_ppm_hat_r.F \
 gad_ppm_hat_x.F \
 gad_ppm_hat_y.F \
 gad_ppm_p3e_r.F \
 gad_ppm_p3e_x.F \
 gad_ppm_p3e_y.F \
 gad_pqm_adv_r.F \
 gad_pqm_adv_x.F \
 gad_pqm_adv_y.F \
 gad_pqm_flx_r.F \
 gad_pqm_flx_x.F \
 gad_pqm_flx_y.F \
 gad_pqm_fun.F \
 gad_pqm_hat_r.F \
 gad_pqm_hat_x.F \
 gad_pqm_hat_y.F \
 gad_pqm_p5e_r.F \
 gad_pqm_p5e_x.F \
 gad_pqm_p5e_y.F \
 gad_read_pickup.F \
 gad_som_adv_r.F \
 gad_som_adv_x.F \
 gad_som_adv_y.F \
 gad_som_advect.F \
 gad_som_exchanges.F \
 gad_som_fill_cs_corner.F \
 gad_som_lim_r.F \
 gad_som_prep_cs_corner.F \
 gad_u3_adv_r.F \
 gad_u3_adv_x.F \
 gad_u3_adv_y.F \
 gad_u3c4_impl_r.F \
 gad_write_pickup.F \
 salt_fill.F \
 mdsio_buffertorl.F \
 mdsio_buffertors.F \
 mdsio_check4file.F \
 mdsio_facef_read.F \
 mdsio_gl.F \
 mdsio_gl_slice.F \
 mdsio_pass_r4torl.F \
 mdsio_pass_r4tors.F \
 mdsio_pass_r8torl.F \
 mdsio_pass_r8tors.F \
 mdsio_rd_rec_rl.F \
 mdsio_rd_rec_rs.F \
 mdsio_read_field.F \
 mdsio_read_meta.F \
 mdsio_read_section.F \
 mdsio_read_tape.F \
 mdsio_read_whalos.F \
 mdsio_readvec_loc.F \
 mdsio_rw_field.F \
 mdsio_rw_slice.F \
 mdsio_seg4torl.F \
 mdsio_seg4tors.F \
 mdsio_seg8torl.F \
 mdsio_seg8tors.F \
 mdsio_segxtorx_2d.F \
 mdsio_wr_metafiles.F \
 mdsio_wr_rec_rl.F \
 mdsio_wr_rec_rs.F \
 mdsio_write_field.F \
 mdsio_write_meta.F \
 mdsio_write_section.F \
 mdsio_write_tape.F \
 mdsio_write_whalos.F \
 mdsio_writelocal.F \
 mdsio_writevec_loc.F \
 mom_calc_3d_strain.F \
 mom_calc_absvort3.F \
 mom_calc_hdiv.F \
 mom_calc_hfacz.F \
 mom_calc_ke.F \
 mom_calc_relvort3.F \
 mom_calc_smag_3d.F \
 mom_calc_strain.F \
 mom_calc_tension.F \
 mom_calc_visc.F \
 mom_diagnostics_init.F \
 mom_hdissip.F \
 mom_init_fixed.F \
 mom_quasihydrostatic.F \
 mom_u_botdrag_coeff.F \
 mom_u_coriolis_nh.F \
 mom_u_implicit_r.F \
 mom_u_metric_nh.F \
 mom_u_rviscflux.F \
 mom_u_sidedrag.F \
 mom_uv_smag_3d.F \
 mom_v_botdrag_coeff.F \
 mom_v_coriolis_nh.F \
 mom_v_implicit_r.F \
 mom_v_metric_nh.F \
 mom_v_rviscflux.F \
 mom_v_sidedrag.F \
 mom_visc_qgl_limit.F \
 mom_visc_qgl_stretch.F \
 mom_w_coriolis_nh.F \
 mom_w_metric_nh.F \
 mom_w_sidedrag.F \
 mom_w_smag_3d.F \
 mom_calc_rtrans.F \
 mom_fluxform.F \
 mom_u_adv_uu.F \
 mom_u_adv_vu.F \
 mom_u_adv_wu.F \
 mom_u_coriolis.F \
 mom_u_del2u.F \
 mom_u_metric_cylinder.F \
 mom_u_metric_sphere.F \
 mom_u_xviscflux.F \
 mom_u_yviscflux.F \
 mom_uv_boundary.F \
 mom_v_adv_uv.F \
 mom_v_adv_vv.F \
 mom_v_adv_wv.F \
 mom_v_coriolis.F \
 mom_v_del2v.F \
 mom_v_metric_cylinder.F \
 mom_v_metric_sphere.F \
 mom_v_xviscflux.F \
 mom_v_yviscflux.F \
 mom_vecinv.F \
 mom_vi_coriolis.F \
 mom_vi_del2uv.F \
 mom_vi_hdissip.F \
 mom_vi_u_coriolis.F \
 mom_vi_u_coriolis_c4.F \
 mom_vi_u_grad_ke.F \
 mom_vi_u_vertshear.F \
 mom_vi_v_coriolis.F \
 mom_vi_v_coriolis_c4.F \
 mom_vi_v_grad_ke.F \
 mom_vi_v_vertshear.F \
 mon_advcfl.F \
 mon_advcflw.F \
 mon_advcflw2.F \
 mon_calc_advcfl.F \
 mon_calc_stats_rl.F \
 mon_calc_stats_rs.F \
 mon_init.F \
 mon_ke.F \
 mon_out.F \
 mon_printstats_rl.F \
 mon_printstats_rs.F \
 mon_set_iounit.F \
 mon_set_pref.F \
 mon_solution.F \
 mon_stats_latbnd_rl.F \
 mon_stats_rl.F \
 mon_stats_rs.F \
 mon_surfcor.F \
 mon_vort3.F \
 mon_writestats_rl.F \
 mon_writestats_rs.F \
 monitor.F \
 monitor_ad.F \
 monitor_g.F \
 mypackage_calc_rhs.F \
 mypackage_check.F \
 mypackage_diagnostics_init.F \
 mypackage_diagnostics_state.F \
 mypackage_init_fixed.F \
 mypackage_init_varia.F \
 mypackage_mnc_init.F \
 mypackage_read_pickup.F \
 mypackage_readparms.F \
 mypackage_tendency_apply.F \
 mypackage_write_pickup.F \
 get_write_global_fld.F \
 read_fld_xy_rl.F \
 read_fld_xy_rs.F \
 read_fld_xyz_rl.F \
 read_fld_xyz_rs.F \
 read_glvec_rl.F \
 read_glvec_rs.F \
 read_mflds.F \
 read_rec.F \
 rw_get_suffix.F \
 set_write_global_fld.F \
 write_fld_3d_rl.F \
 write_fld_3d_rs.F \
 write_fld_s3d_rl.F \
 write_fld_s3d_rs.F \
 write_fld_xy_rl.F \
 write_fld_xy_rs.F \
 write_fld_xyz_rl.F \
 write_fld_xyz_rs.F \
 write_glvec_rl.F \
 write_glvec_rs.F \
 write_local_rl.F \
 write_local_rs.F \
 write_rec.F \
 all_proc_die.F \
 bar2.F \
 bar_check.F \
 barrier.F \
 check_threads.F \
 comm_stats.F \
 cumulsum_z_tile.F \
 diff_phase_multiple.F \
 different_multiple.F \
 eeboot.F \
 eeboot_minimal.F \
 eedata_example.F \
 eedie.F \
 eeintro_msg.F \
 eeset_parms.F \
 eewrite_eeenv.F \
 exch0_r4.F \
 exch0_r8.F \
 exch0_rl.F \
 exch0_rs.F \
 exch1_bg_r4_cube.F \
 exch1_bg_r8_cube.F \
 exch1_bg_rl_cube.F \
 exch1_bg_rs_cube.F \
 exch1_r4.F \
 exch1_r4_cube.F \
 exch1_r8.F \
 exch1_r8_cube.F \
 exch1_rl.F \
 exch1_rl_ad.F \
 exch1_rl_cube.F \
 exch1_rl_cube_ad.F \
 exch1_rs.F \
 exch1_rs_ad.F \
 exch1_rs_cube.F \
 exch1_rs_cube_ad.F \
 exch1_uv_r4_cube.F \
 exch1_uv_r8_cube.F \
 exch1_uv_rl_cube.F \
 exch1_uv_rs_cube.F \
 exch1_z_r4_cube.F \
 exch1_z_r8_cube.F \
 exch1_z_rl_cube.F \
 exch1_z_rs_cube.F \
 exch_3d_r4.F \
 exch_3d_r8.F \
 exch_3d_rl.F \
 exch_3d_rs.F \
 exch_cycle_ebl.F \
 exch_init.F \
 exch_r4_recv_get_x.F \
 exch_r4_recv_get_y.F \
 exch_r4_send_put_x.F \
 exch_r4_send_put_y.F \
 exch_r8_recv_get_x.F \
 exch_r8_recv_get_y.F \
 exch_r8_send_put_x.F \
 exch_r8_send_put_y.F \
 exch_rl_recv_get_x.F \
 exch_rl_recv_get_y.F \
 exch_rl_send_put_x.F \
 exch_rl_send_put_y.F \
 exch_rs_recv_get_x.F \
 exch_rs_recv_get_y.F \
 exch_rs_send_put_x.F \
 exch_rs_send_put_y.F \
 exch_s3d_r4.F \
 exch_s3d_r8.F \
 exch_s3d_rl.F \
 exch_s3d_rs.F \
 exch_sm_3d_r4.F \
 exch_sm_3d_r8.F \
 exch_sm_3d_rl.F \
 exch_sm_3d_rs.F \
 exch_tap_b.F \
 exch_tap_d.F \
 exch_uv_3d_r4.F \
 exch_uv_3d_r8.F \
 exch_uv_3d_rl.F \
 exch_uv_3d_rs.F \
 exch_uv_agrid_3d_r4.F \
 exch_uv_agrid_3d_r8.F \
 exch_uv_agrid_3d_rl.F \
 exch_uv_agrid_3d_rs.F \
 exch_uv_bgrid_3d_r4.F \
 exch_uv_bgrid_3d_r8.F \
 exch_uv_bgrid_3d_rl.F \
 exch_uv_bgrid_3d_rs.F \
 exch_uv_dgrid_3d_r4.F \
 exch_uv_dgrid_3d_r8.F \
 exch_uv_dgrid_3d_rl.F \
 exch_uv_dgrid_3d_rs.F \
 exch_uv_xy_r4.F \
 exch_uv_xy_r8.F \
 exch_uv_xy_rl.F \
 exch_uv_xy_rs.F \
 exch_uv_xyz_r4.F \
 exch_uv_xyz_r8.F \
 exch_uv_xyz_rl.F \
 exch_uv_xyz_rs.F \
 exch_xy_r4.F \
 exch_xy_r8.F \
 exch_xy_rl.F \
 exch_xy_rs.F \
 exch_xyz_r4.F \
 exch_xyz_r8.F \
 exch_xyz_rl.F \
 exch_xyz_rs.F \
 exch_z_3d_r4.F \
 exch_z_3d_r8.F \
 exch_z_3d_rl.F \
 exch_z_3d_rs.F \
 fill_cs_corner_ag_rl.F \
 fill_cs_corner_tr_rl.F \
 fill_cs_corner_uv_rl.F \
 fill_cs_corner_uv_rs.F \
 fool_the_compiler.F \
 gather_2d_r4.F \
 gather_2d_r8.F \
 gather_2d_wh_r4.F \
 gather_2d_wh_r8.F \
 gather_vec_r4.F \
 gather_vec_r8.F \
 gather_xz.F \
 gather_yz.F \
 get_periodic_interval.F \
 global_max.F \
 global_sum.F \
 global_sum_singlecpu.F \
 global_sum_tile.F \
 global_sum_vector.F \
 ini_communication_patterns.F \
 ini_procs.F \
 ini_threading_environment.F \
 lef_zero.F \
 main.F \
 master_cpu_io.F \
 master_cpu_thread.F \
 mds_byteswapi4.F \
 mds_byteswapr4.F \
 mds_byteswapr8.F \
 mds_flush.F \
 mds_reclen.F \
 mdsfindunit.F \
 memsync.F \
 nml_change_syntax.F \
 nml_set_terminator.F \
 open_copy_data_file.F \
 print.F \
 reset_halo.F \
 scatter_2d_r4.F \
 scatter_2d_r8.F \
 scatter_2d_wh_r4.F \
 scatter_2d_wh_r8.F \
 scatter_vec_r4.F \
 scatter_vec_r8.F \
 scatter_xz.F \
 scatter_yz.F \
 stop_if_error.F \
 timers.F \
 utils.F \
 write_utils.F \
 adams_bashforth2.F \
 adams_bashforth3.F \
 add_walls2masks.F \
 apply_forcing.F \
 calc_3d_diffusivity.F \
 calc_adv_flow.F \
 calc_div_ghat.F \
 calc_eddy_stress.F \
 calc_grad_phi_fv.F \
 calc_grad_phi_hyd.F \
 calc_grad_phi_surf.F \
 calc_grid_angles.F \
 calc_gw.F \
 calc_ivdc.F \
 calc_oce_mxlayer.F \
 calc_phi_hyd.F \
 calc_r_star.F \
 calc_surf_dr.F \
 calc_viscosity.F \
 calc_wsurf_tr.F \
 cg2d.F \
 cg2d_ex0.F \
 cg2d_nsa.F \
 cg2d_sr.F \
 cg3d.F \
 cg3d_ex0.F \
 check_pickup.F \
 config_check.F \
 config_summary.F \
 convective_adjustment.F \
 convective_adjustment_ini.F \
 convective_weights.F \
 convectively_mixtracer.F \
 correction_step.F \
 cycle_ab_tracer.F \
 cycle_tracer.F \
 diags_oceanic_surf_flux.F \
 diags_phi_hyd.F \
 diags_phi_rlow.F \
 diags_rho.F \
 diags_sound_speed.F \
 do_atmospheric_phys.F \
 do_fields_blocking_exchanges.F \
 do_oceanic_phys.F \
 do_stagger_fields_exchanges.F \
 do_statevars_diags.F \
 do_the_model_io.F \
 do_write_pickup.F \
 dynamics.F \
 external_fields_load.F \
 external_forcing.F \
 external_forcing_surf.F \
 find_alpha.F \
 find_hyd_press_1d.F \
 find_rho.F \
 forcing_surf_relax.F \
 forward_step.F \
 freesurf_rescale_g.F \
 freeze_surface.F \
 grad_sigma.F \
 gsw_teos10.F \
 impldiff.F \
 ini_cartesian_grid.F \
 ini_cg2d.F \
 ini_cg3d.F \
 ini_cori.F \
 ini_curvilinear_grid.F \
 ini_cylinder_grid.F \
 ini_depths.F \
 ini_dynvars.F \
 ini_eos.F \
 ini_ffields.F \
 ini_fields.F \
 ini_forcing.F \
 ini_global_domain.F \
 ini_grid.F \
 ini_linear_phisurf.F \
 ini_local_grid.F \
 ini_masks_etc.F \
 ini_mixing.F \
 ini_mnc_vars.F \
 ini_model_io.F \
 ini_nh_fields.F \
 ini_nh_vars.F \
 ini_nlfs_vars.F \
 ini_p_ground.F \
 ini_parms.F \
 ini_pressure.F \
 ini_psurf.F \
 ini_salt.F \
 ini_sigma_hfac.F \
 ini_spherical_polar_grid.F \
 ini_theta.F \
 ini_vel.F \
 ini_vertical_grid.F \
 initialise_fixed.F \
 initialise_varia.F \
 integr_continuity.F \
 integrate_for_w.F \
 load_fields_driver.F \
 load_grid_spacing.F \
 load_ref_files.F \
 main_do_loop.F \
 momentum_correction_step.F \
 packages_boot.F \
 packages_check.F \
 packages_error_msg.F \
 packages_init_fixed.F \
 packages_init_variables.F \
 packages_print_msg.F \
 packages_readparms.F \
 packages_unused_msg.F \
 packages_write_pickup.F \
 plot_field.F \
 port_rand.F \
 post_cg3d.F \
 pre_cg3d.F \
 pressure_for_eos.F \
 read_pickup.F \
 remove_mean.F \
 reset_nlfs_vars.F \
 rotate_spherical_polar_grid.F \
 rotate_uv2en.F \
 salt_integrate.F \
 seawater.F \
 set_defaults.F \
 set_grid_factors.F \
 set_parms.F \
 set_ref_state.F \
 solve_for_pressure.F \
 solve_pentadiagonal.F \
 solve_tridiagonal.F \
 solve_uv_tridiago.F \
 state_summary.F \
 swfrac.F \
 taueddy_init_varia.F \
 taueddy_tendency_apply.F \
 temp_integrate.F \
 the_main_loop.F \
 the_model_main.F \
 thermodynamics.F \
 timestep.F \
 timestep_tracer.F \
 timestep_wvel.F \
 tracers_correction_step.F \
 tracers_iigw_correction.F \
 turnoff_model_io.F \
 update_cg2d.F \
 update_etah.F \
 update_etaws.F \
 update_masks_etc.F \
 update_r_star.F \
 update_sigma.F \
 update_surf_dr.F \
 write_grid.F \
 write_pickup.F \
 write_state.F
F90_SRC_FILES = 
C_SRC_FILES =  \
 gsl_ieee_env.c \
 ptwrapper.c \
 setdir.c \
 setrlstk.c \
 sigreg.c \
 tim.c \
 timer_stats.c
H_SRC_FILES =  \
 DIAGNOSTICS_SIZE.h \
 SIZE.h \
 DEBUG_OPTIONS.h \
 DIAGNOSTICS.h \
 DIAGNOSTICS_CALC.h \
 DIAGNOSTICS_P2SHARE.h \
 DIAGSTATS_REGIONS.h \
 DIAG_OPTIONS.h \
 GAD.h \
 GAD_FLUX_LIMITER.h \
 GAD_OPTIONS.h \
 GAD_SOM_VARS.h \
 gad_ad_check_lev1_dir.h \
 gad_ad_check_lev2_dir.h \
 gad_ad_check_lev3_dir.h \
 gad_ad_check_lev4_dir.h \
 MDSIO_BUFF_3D.h \
 MDSIO_BUFF_WH.h \
 MDSIO_OPTIONS.h \
 MOM_COMMON_OPTIONS.h \
 MOM_VISC.h \
 MOM_FLUXFORM.h \
 MOM_FLUXFORM_OPTIONS.h \
 MOM_VECINV_OPTIONS.h \
 MONITOR.h \
 MONITOR_OPTIONS.h \
 MYPACKAGE.h \
 MYPACKAGE_OPTIONS.h \
 RW_MFLDS.h \
 RW_OPTIONS.h \
 mpi-ext.h \
 mpi.h \
 mpi_portable_platform.h \
 mpif-c-constants-decl.h \
 mpif-config.h \
 mpif-constants.h \
 mpif-ext.h \
 mpif-externals.h \
 mpif-handles.h \
 mpif-io-constants.h \
 mpif-io-handles.h \
 mpif-sentinels.h \
 mpif-sizeof.h \
 mpif.h \
 pshmem.h \
 pshmemx.h \
 shmem-compat.h \
 shmem.h \
 shmemx.h \
 BAR2.h \
 BARRIER.h \
 CPP_EEMACROS.h \
 CPP_EEOPTIONS.h \
 CUMULSUM.h \
 DEF_IN_MAKEFILE.h \
 EEBUFF_SCPU.h \
 EEPARAMS.h \
 EESUPPORT.h \
 EXCH.h \
 GLOBAL_MAX.h \
 GLOBAL_SUM.h \
 MAIN_PDIRECTIVES1.h \
 MAIN_PDIRECTIVES2.h \
 MPI_INFO.h \
 SIGREG.h \
 CG2D.h \
 CG3D.h \
 CPP_OPTIONS.h \
 DXC_MACROS.h \
 DXF_MACROS.h \
 DXG_MACROS.h \
 DXV_MACROS.h \
 DYC_MACROS.h \
 DYF_MACROS.h \
 DYG_MACROS.h \
 DYNVARS.h \
 DYU_MACROS.h \
 EOS.h \
 FCORI_MACROS.h \
 FFIELDS.h \
 GRID.h \
 GRID_MACROS.h \
 HFACC_MACROS.h \
 HFACS_MACROS.h \
 HFACW_MACROS.h \
 MASKS_MACROS.h \
 MASKW_MACROS.h \
 NH_VARS.h \
 PARAMS.h \
 RAS_MACROS.h \
 RAW_MACROS.h \
 RA_MACROS.h \
 RECIP_DXC_MACROS.h \
 RECIP_DXF_MACROS.h \
 RECIP_DXG_MACROS.h \
 RECIP_DXV_MACROS.h \
 RECIP_DYC_MACROS.h \
 RECIP_DYF_MACROS.h \
 RECIP_DYG_MACROS.h \
 RECIP_DYU_MACROS.h \
 RECIP_HFACC_MACROS.h \
 RECIP_HFACS_MACROS.h \
 RECIP_HFACW_MACROS.h \
 RESTART.h \
 SET_GRID.h \
 SURFACE.h \
 TANPHIATU_MACROS.h \
 TANPHIATV_MACROS.h \
 XC_MACROS.h \
 YC_MACROS.h
AD_FLOW_FILES =  \
 debug_ad.flow \
 diagnostics_ad.flow \
 mdsio_ad.flow \
 monitor_ad.flow \
 mypackage_ad.flow \
 rw_ad.flow \
 eesupp.flow \
 exch1_ad.flow

F77_PP_SRC_FILES = $(F77_SRC_FILES:.F=.f)
F90_PP_SRC_FILES = $(F90_SRC_FILES:.F90=.f90)
OBJFILES= $(F90_SRC_FILES:.F90=.o) $(F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o)

.SUFFIXES:
.SUFFIXES: .o .f .p .F .c .ff90 .f90 .FF90 .F90 .flow

all: fwd_exe_target
fwd_exe_target:
	@echo Update AD_CONFIG.h and make $(EXECUTABLE)
	@$(TOOLSDIR)/convert_cpp_cmd2defines "Forward version" -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -UALLOW_TANGENTLINEAR_RUN > ad_config.template
	@cmp ad_config.template AD_CONFIG.h || cat ad_config.template > AD_CONFIG.h
	@-rm -f ad_config.template
	$(MAKE) -f $(MAKEFILE) $(EXECUTABLE)

$(EXECUTABLE): $(SPECIAL_FILES) $(F90_SRC_FILES) $(F77_SRC_FILES) $(C_SRC_FILES) $(H_SRC_FILES) $(OBJFILES) $(EMBEDDED_FILES)
	@echo Creating $@ ...
	$(LINK) -o $@ $(FFLAGS) $(FOPTIM) $(OBJFILES) $(LIBS)

depend:
	@$(MAKE) -f $(MAKEFILE) links
	$(MAKEDEPEND) -f $(MAKEFILE) -o .f $(DEFINES) $(INCLUDES) $(CPPINCLUDES) $(F77_SRC_FILES)
	$(MAKEDEPEND) -f $(MAKEFILE) -a -o .f90 $(DEFINES) $(INCLUDES) $(CPPINCLUDES) $(F90_SRC_FILES)
	$(TOOLSDIR)/f90mkdepend -fs f -fs90 f90 >> $(MAKEFILE)
	-rm -f makedepend.out

lib: libmitgcmuv.a

libmitgcmuv.a: $(SPECIAL_FILES) $(OBJFILES)
	ar rcv libmitgcmuv.a $(OBJFILES)
	ar d   libmitgcmuv.a main.o

obj: $(OBJFILES)

links: $(F90_SRC_FILES) $(F77_SRC_FILES) $(C_SRC_FILES) $(H_SRC_FILES) $(SPECIAL_FILES)

small_f: $(F90_PP_SRC_FILES) $(F77_PP_SRC_FILES)

# remove most of the files that "make" generates
clean:
	-rm -rf *.p *.f90 *.mod  work.{pc,pcl} *.template
	-rm -rf *.o
	-rm -rf *.f
	-rm -rf *.ff90 $(AD_CLEAN) ad_input*

# remove most of the files that "make" and "make depend" generate
Clean:
	@$(MAKE) -f $(MAKEFILE) clean
	@$(MAKE) -f $(MAKEFILE) cleanlinks
	-rm -rf *.flowdir
	-rm -f $(SPECIAL_FILES) f90mkdepend.log Makefile.old
	-rm -f taf_command taf_output taf_ad.log taf_ad_flow.log taf_ftl.log

# remove also the executable, files that "genmake2" generates (except Makefile)
#         and output from a run (plus log files from testreport)
CLEAN:
	@$(MAKE) -f $(MAKEFILE) Clean
	-rm -f $(EXECUTABLE) $(EXE_ADJ) $(EXE_TLM) *.bak
	-rm -f genmake.log genmake_state genmake_*optfile
	-rm -f SIZE.h.mpi genmake.tr_log make.tr_log
	-rm -rf mpi_headers mnc_test_*
	-find $(EXEDIR) -name "*.meta" -exec rm {} \;
	-find $(EXEDIR) -name "*.data" -exec rm {} \;
	-find $(EXEDIR) -name "fort.*" -exec rm {} \;
	-rm -f *.txt STDOUT.* STDERR.* *diagnostics.log *.[0-9][0-9][0-9][0-9].log
	-rm -f datetime costfinal divided.ctrl snapshot* output_adm.txt.diva_*
	-rm -f *_MIT_CE_000.opt0000 costfunction*0000
	-rm -f oad_cp.[0-9][0-9][0-9].?????

makefile:
	$(TOOLSDIR)/genmake2  "-rootdir=/home/hmh85/scratch/MITgcm" "-mpi" "-mods=/home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code" "-of=/home/hmh85/scratch/MITgcm/tools/build_options/linux_amd64_gfortran"
cleanlinks:
	-find . -type l -exec rm {} \;

# Special targets (SPECIAL_FILES) which are created by make
PACKAGES_CONFIG.h:
	@echo Creating $@ ...
	@$(TOOLSDIR)/convert_cpp_cmd2defines -bPACKAGES_CONFIG_H "Disabled packages:" $(DISABLED_PACKAGES) " " "Enabled packages:" $(ENABLED_PACKAGES) > $@
AD_CONFIG.h:
	@echo Creating $@ ...
	@$(TOOLSDIR)/convert_cpp_cmd2defines "Default version" -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -UALLOW_TANGENTLINEAR_RUN > $@
FC_NAMEMANGLE.h:
	@echo Creating $@ ...
	echo "#define FC_NAMEMANGLE(X)  X ## _" > $@

BUILD_INFO.h:
	@echo Creating $@ ...
	-echo "#define THISVER 'checkpoint69m'" > $@
	-echo "#define THISUSER 'hmh85'" >> $@
	-echo "#define THISDATE 'Thu Jun  4 17:46:47 EEST 2026'" >> $@
	-echo "#define THISHOST 'anode07'" >> $@

# The normal chain of rules is (  .F - .f - .o  )

## This nullifies any default implicit rules concerning these two file types:
## %.o : %.F

# C preprocessing and replacing the _d in constants:
CPPCMD = cat $< |  cpp -traditional -P $(DEFINES) $(INCLUDES) | $(TOOLSDIR)/set64bitConst.sh

.F.f:
	$(CPPCMD)  > $@
.f.o:
	$(FC) $(FFLAGS) $(FOPTIM) -c $<
.F.o:
	$(FC) $(FFLAGS) $(FOPTIM) -c $<
.F90.f90:
	$(CPPCMD)  > $@
.FF90.ff90:
	$(CPPCMD)  > $@
.f90.o:
	$(F90C) $(F90FLAGS) $(F90OPTIM) -c $<
.ff90.o:
	cp $< $(basename  $<).f90
	$(F90C) $(F90FLAGS) $(F90OPTIM) -c $(F90FIXEDFORMAT) $(basename  $<).f90
.c.o:
	$(CC) $(CFLAGS) $(DEFINES) $(INCLUDES) -c $<

# Special exceptions that use the ( .F - .p - .f - .o ) rule-chain
.F.p:
	$(CPPCMD) > $@
.p.f:
	$(KPP) $(KFLAGS1)$@ $(KFLAGS2) $<


TAF           = staf
TAF_EXTRA     = 
TAF_F77_FLAGS = 
TAF_F90_FLAGS = 
LOCAL_MPI_HEADERS = 

AD_TAF_FLAGS  = -f77 -i4 -r4 -intrinsic system,flush -reverse -l taf_ad.log -toplevel 'the_main_loop' -input 'xx_obcsn_dummy xx_obcss_dummy xx_obcsw_dummy xx_obcse_dummy xx_genarr2d_dummy xx_genarr3d_dummy xx_gentim2d_dummy' -output 'fc'
FTL_TAF_FLAGS = -f77 -i4 -r4 -intrinsic system,flush -forward -l taf_ftl.log -toplevel 'the_main_loop' -input 'xx_obcsn_dummy xx_obcss_dummy xx_obcsw_dummy xx_obcse_dummy xx_genarr2d_dummy xx_genarr3d_dummy xx_gentim2d_dummy' -output 'fc'
SVD_TAF_FLAGS = -f77 -i4 -r4 -intrinsic system,flush -reverse -forward -pure -l taf_svd.log -toplevel 'the_main_loop' -input 'xx_obcsn_dummy xx_obcss_dummy xx_obcsw_dummy xx_obcse_dummy xx_genarr2d_dummy xx_genarr3d_dummy xx_gentim2d_dummy' -output 'fc'

AD_FILES =  \
 diagnostics_is_on.f \
 gad_advection.f \
 gad_biharm_r.f \
 gad_biharm_x.f \
 gad_biharm_y.f \
 gad_c2_adv_r.f \
 gad_c2_impl_r.f \
 gad_c2_adv_x.f \
 gad_c2_adv_y.f \
 gad_c4_adv_r.f \
 gad_c4_adv_x.f \
 gad_c4_adv_y.f \
 gad_calc_rhs.f \
 gad_del2.f \
 gad_diff_r.f \
 gad_diff_x.f \
 gad_diff_y.f \
 gad_dst2u1_adv_r.f \
 gad_dst2u1_impl_r.f \
 gad_dst2u1_adv_x.f \
 gad_dst2u1_adv_y.f \
 gad_dst3_adv_x.f \
 gad_dst3_adv_y.f \
 gad_dst3_adv_r.f \
 gad_dst3fl_adv_x.f \
 gad_dst3fl_adv_y.f \
 gad_dst3fl_adv_r.f \
 gad_dst3fl_impl_r.f \
 gad_fluxlimit_adv_r.f \
 gad_fluxlimit_impl_r.f \
 gad_fluxlimit_adv_x.f \
 gad_fluxlimit_adv_y.f \
 gad_grad_x.f \
 gad_grad_y.f \
 gad_implicit_r.f \
 gad_init_varia.f \
 gad_os7mp_adv_x.f \
 gad_os7mp_adv_y.f \
 gad_os7mp_adv_r.f \
 gad_read_pickup.f \
 gad_u3_adv_r.f \
 gad_u3c4_impl_r.f \
 gad_u3_adv_x.f \
 gad_u3_adv_y.f \
 gad_som_advect.f \
 gad_som_adv_r.f \
 gad_som_adv_x.f \
 gad_som_adv_y.f \
 gad_som_fill_cs_corner.f \
 gad_som_lim_r.f \
 gad_som_prep_cs_corner.f \
 gad_som_exchanges.f \
 gad_exch_som.f \
 mom_calc_absvort3.f \
 mom_calc_hdiv.f \
 mom_calc_hfacz.f \
 mom_calc_ke.f \
 mom_calc_relvort3.f \
 mom_calc_strain.f \
 mom_calc_tension.f \
 mom_calc_visc.f \
 mom_hdissip.f \
 mom_quasihydrostatic.f \
 mom_u_botdrag_coeff.f \
 mom_u_coriolis_nh.f \
 mom_u_metric_nh.f \
 mom_u_rviscflux.f \
 mom_u_sidedrag.f \
 mom_v_botdrag_coeff.f \
 mom_v_coriolis_nh.f \
 mom_v_metric_nh.f \
 mom_v_rviscflux.f \
 mom_v_sidedrag.f \
 mom_visc_qgl_limit.f \
 mom_visc_qgl_stretch.f \
 mom_calc_rtrans.f \
 mom_fluxform.f \
 mom_u_adv_uu.f \
 mom_u_adv_vu.f \
 mom_u_adv_wu.f \
 mom_u_coriolis.f \
 mom_u_del2u.f \
 mom_u_metric_cylinder.f \
 mom_u_metric_sphere.f \
 mom_u_xviscflux.f \
 mom_u_yviscflux.f \
 mom_v_adv_uv.f \
 mom_v_adv_vv.f \
 mom_v_adv_wv.f \
 mom_v_coriolis.f \
 mom_v_del2v.f \
 mom_v_metric_cylinder.f \
 mom_v_metric_sphere.f \
 mom_v_xviscflux.f \
 mom_v_yviscflux.f \
 mom_vecinv.f \
 mom_vi_coriolis.f \
 mom_vi_del2uv.f \
 mom_vi_hdissip.f \
 mom_vi_u_coriolis.f \
 mom_vi_u_coriolis_c4.f \
 mom_vi_u_grad_ke.f \
 mom_vi_u_vertshear.f \
 mom_vi_v_coriolis.f \
 mom_vi_v_coriolis_c4.f \
 mom_vi_v_grad_ke.f \
 mom_vi_v_vertshear.f \
 mypackage_init_varia.f \
 mypackage_tendency_apply.f \
 mypackage_calc_rhs.f \
 read_fld_xyz_rs.f \
 read_fld_xyz_rl.f \
 read_fld_xy_rs.f \
 read_fld_xy_rl.f \
 read_rec.f \
 different_multiple.f \
 fill_cs_corner_ag_rl.f \
 fill_cs_corner_tr_rl.f \
 fill_cs_corner_uv_rs.f \
 fill_cs_corner_uv_rl.f \
 get_periodic_interval.f \
 lef_zero.f \
 exch_sm_3d_rl.f \
 exch_3d_rs.f \
 exch_3d_rl.f \
 exch_xy_rl.f \
 exch_xy_rs.f \
 exch_xyz_rl.f \
 exch_xyz_rs.f \
 exch_uv_3d_rs.f \
 exch_uv_3d_rl.f \
 exch_uv_xy_rs.f \
 exch_uv_xy_rl.f \
 exch_uv_xyz_rs.f \
 exch_uv_xyz_rl.f \
 exch_uv_agrid_3d_rs.f \
 exch_uv_agrid_3d_rl.f \
 exch_uv_dgrid_3d_rs.f \
 exch_uv_dgrid_3d_rl.f \
 adams_bashforth2.f \
 adams_bashforth3.f \
 apply_forcing.f \
 calc_3d_diffusivity.f \
 calc_adv_flow.f \
 calc_div_ghat.f \
 calc_grad_phi_fv.f \
 calc_grad_phi_hyd.f \
 calc_grad_phi_surf.f \
 calc_ivdc.f \
 calc_phi_hyd.f \
 calc_r_star.f \
 calc_surf_dr.f \
 calc_viscosity.f \
 calc_wsurf_tr.f \
 cg2d_nsa.f \
 convective_adjustment.f \
 convective_adjustment_ini.f \
 convective_weights.f \
 convectively_mixtracer.f \
 correction_step.f \
 cycle_ab_tracer.f \
 cycle_tracer.f \
 diags_phi_hyd.f \
 diags_phi_rlow.f \
 diags_rho.f \
 diags_sound_speed.f \
 do_atmospheric_phys.f \
 do_fields_blocking_exchanges.f \
 do_oceanic_phys.f \
 do_stagger_fields_exchanges.f \
 dynamics.f \
 external_fields_load.f \
 external_forcing_surf.f \
 forcing_surf_relax.f \
 find_alpha.f \
 find_rho.f \
 forward_step.f \
 freesurf_rescale_g.f \
 freeze_surface.f \
 grad_sigma.f \
 impldiff.f \
 ini_cg2d.f \
 ini_depths.f \
 ini_ffields.f \
 ini_fields.f \
 ini_forcing.f \
 ini_mixing.f \
 ini_nh_fields.f \
 ini_nlfs_vars.f \
 ini_p_ground.f \
 ini_pressure.f \
 ini_psurf.f \
 ini_salt.f \
 ini_theta.f \
 ini_dynvars.f \
 ini_vel.f \
 initialise_varia.f \
 integrate_for_w.f \
 integr_continuity.f \
 load_fields_driver.f \
 load_grid_spacing.f \
 main_do_loop.f \
 momentum_correction_step.f \
 packages_init_variables.f \
 pressure_for_eos.f \
 read_pickup.f \
 remove_mean.f \
 reset_nlfs_vars.f \
 rotate_uv2en.f \
 salt_integrate.f \
 temp_integrate.f \
 seawater.f \
 solve_for_pressure.f \
 state_summary.f \
 solve_pentadiagonal.f \
 solve_tridiagonal.f \
 taueddy_init_varia.f \
 taueddy_tendency_apply.f \
 timestep.f \
 the_main_loop.f \
 thermodynamics.f \
 timestep_tracer.f \
 tracers_correction_step.f \
 update_cg2d.f \
 update_etah.f \
 update_etaws.f \
 update_masks_etc.f \
 update_r_star.f \
 update_sigma.f \
 update_surf_dr.f \
 swfrac.f
AD_F90FILES = 
NON_AD_F77_SRC_FILES =  \
 chksum_tiled.F \
 debug_call.F \
 debug_cs_corner_uv.F \
 debug_enter.F \
 debug_fld_stats_rl.F \
 debug_fld_stats_rs.F \
 debug_leave.F \
 debug_msg.F \
 debug_stats_rl.F \
 debug_stats_rs.F \
 fill_in_corners_rl.F \
 write_fullarray_rl.F \
 write_fullarray_rs.F \
 diag_calc_psivel.F \
 diag_cg2d.F \
 diag_vegtile_fill.F \
 diagnostics_addtolist.F \
 diagnostics_calc_phivel.F \
 diagnostics_check.F \
 diagnostics_clear.F \
 diagnostics_fill.F \
 diagnostics_fill_field.F \
 diagnostics_fill_rs.F \
 diagnostics_fill_state.F \
 diagnostics_fract_fill.F \
 diagnostics_ini_io.F \
 diagnostics_init_early.F \
 diagnostics_init_fixed.F \
 diagnostics_init_varia.F \
 diagnostics_interp_p2p.F \
 diagnostics_interp_vert.F \
 diagnostics_list_check.F \
 diagnostics_main_init.F \
 diagnostics_mnc_out.F \
 diagnostics_out.F \
 diagnostics_read_pickup.F \
 diagnostics_readparms.F \
 diagnostics_scale_fill.F \
 diagnostics_scale_fill_rs.F \
 diagnostics_set_calc.F \
 diagnostics_set_levels.F \
 diagnostics_set_pointers.F \
 diagnostics_setdiag.F \
 diagnostics_status_error.F \
 diagnostics_sum_levels.F \
 diagnostics_summary.F \
 diagnostics_switch_onoff.F \
 diagnostics_utils.F \
 diagnostics_write.F \
 diagnostics_write_adj.F \
 diagnostics_write_pickup.F \
 diagstats_ascii_out.F \
 diagstats_calc.F \
 diagstats_clear.F \
 diagstats_close_io.F \
 diagstats_fill.F \
 diagstats_global.F \
 diagstats_ini_io.F \
 diagstats_local.F \
 diagstats_mnc_out.F \
 diagstats_others_calc.F \
 diagstats_output.F \
 diagstats_set_pointers.F \
 diagstats_set_regions.F \
 diagstats_setdiag.F \
 gad_advscheme.F \
 gad_check.F \
 gad_diagnostics_init.F \
 gad_diagnostics_state.F \
 gad_init_fixed.F \
 gad_osc_hat_r.F \
 gad_osc_hat_x.F \
 gad_osc_hat_y.F \
 gad_osc_mul_r.F \
 gad_osc_mul_x.F \
 gad_osc_mul_y.F \
 gad_plm_fun.F \
 gad_ppm_adv_r.F \
 gad_ppm_adv_x.F \
 gad_ppm_adv_y.F \
 gad_ppm_flx_r.F \
 gad_ppm_flx_x.F \
 gad_ppm_flx_y.F \
 gad_ppm_fun.F \
 gad_ppm_hat_r.F \
 gad_ppm_hat_x.F \
 gad_ppm_hat_y.F \
 gad_ppm_p3e_r.F \
 gad_ppm_p3e_x.F \
 gad_ppm_p3e_y.F \
 gad_pqm_adv_r.F \
 gad_pqm_adv_x.F \
 gad_pqm_adv_y.F \
 gad_pqm_flx_r.F \
 gad_pqm_flx_x.F \
 gad_pqm_flx_y.F \
 gad_pqm_fun.F \
 gad_pqm_hat_r.F \
 gad_pqm_hat_x.F \
 gad_pqm_hat_y.F \
 gad_pqm_p5e_r.F \
 gad_pqm_p5e_x.F \
 gad_pqm_p5e_y.F \
 gad_write_pickup.F \
 salt_fill.F \
 mdsio_buffertorl.F \
 mdsio_buffertors.F \
 mdsio_check4file.F \
 mdsio_facef_read.F \
 mdsio_gl.F \
 mdsio_gl_slice.F \
 mdsio_pass_r4torl.F \
 mdsio_pass_r4tors.F \
 mdsio_pass_r8torl.F \
 mdsio_pass_r8tors.F \
 mdsio_rd_rec_rl.F \
 mdsio_rd_rec_rs.F \
 mdsio_read_field.F \
 mdsio_read_meta.F \
 mdsio_read_section.F \
 mdsio_read_tape.F \
 mdsio_read_whalos.F \
 mdsio_readvec_loc.F \
 mdsio_rw_field.F \
 mdsio_rw_slice.F \
 mdsio_seg4torl.F \
 mdsio_seg4tors.F \
 mdsio_seg8torl.F \
 mdsio_seg8tors.F \
 mdsio_segxtorx_2d.F \
 mdsio_wr_metafiles.F \
 mdsio_wr_rec_rl.F \
 mdsio_wr_rec_rs.F \
 mdsio_write_field.F \
 mdsio_write_meta.F \
 mdsio_write_section.F \
 mdsio_write_tape.F \
 mdsio_write_whalos.F \
 mdsio_writelocal.F \
 mdsio_writevec_loc.F \
 mom_calc_3d_strain.F \
 mom_calc_smag_3d.F \
 mom_diagnostics_init.F \
 mom_init_fixed.F \
 mom_u_implicit_r.F \
 mom_uv_smag_3d.F \
 mom_v_implicit_r.F \
 mom_w_coriolis_nh.F \
 mom_w_metric_nh.F \
 mom_w_sidedrag.F \
 mom_w_smag_3d.F \
 mom_uv_boundary.F \
 mon_advcfl.F \
 mon_advcflw.F \
 mon_advcflw2.F \
 mon_calc_advcfl.F \
 mon_calc_stats_rl.F \
 mon_calc_stats_rs.F \
 mon_init.F \
 mon_ke.F \
 mon_out.F \
 mon_printstats_rl.F \
 mon_printstats_rs.F \
 mon_set_iounit.F \
 mon_set_pref.F \
 mon_solution.F \
 mon_stats_latbnd_rl.F \
 mon_stats_rl.F \
 mon_stats_rs.F \
 mon_surfcor.F \
 mon_vort3.F \
 mon_writestats_rl.F \
 mon_writestats_rs.F \
 monitor.F \
 monitor_ad.F \
 monitor_g.F \
 mypackage_check.F \
 mypackage_diagnostics_init.F \
 mypackage_diagnostics_state.F \
 mypackage_init_fixed.F \
 mypackage_mnc_init.F \
 mypackage_read_pickup.F \
 mypackage_readparms.F \
 mypackage_write_pickup.F \
 get_write_global_fld.F \
 read_glvec_rl.F \
 read_glvec_rs.F \
 read_mflds.F \
 rw_get_suffix.F \
 set_write_global_fld.F \
 write_fld_3d_rl.F \
 write_fld_3d_rs.F \
 write_fld_s3d_rl.F \
 write_fld_s3d_rs.F \
 write_fld_xy_rl.F \
 write_fld_xy_rs.F \
 write_fld_xyz_rl.F \
 write_fld_xyz_rs.F \
 write_glvec_rl.F \
 write_glvec_rs.F \
 write_local_rl.F \
 write_local_rs.F \
 write_rec.F \
 all_proc_die.F \
 bar2.F \
 bar_check.F \
 barrier.F \
 check_threads.F \
 comm_stats.F \
 cumulsum_z_tile.F \
 diff_phase_multiple.F \
 eeboot.F \
 eeboot_minimal.F \
 eedata_example.F \
 eedie.F \
 eeintro_msg.F \
 eeset_parms.F \
 eewrite_eeenv.F \
 exch0_r4.F \
 exch0_r8.F \
 exch0_rl.F \
 exch0_rs.F \
 exch1_bg_r4_cube.F \
 exch1_bg_r8_cube.F \
 exch1_bg_rl_cube.F \
 exch1_bg_rs_cube.F \
 exch1_r4.F \
 exch1_r4_cube.F \
 exch1_r8.F \
 exch1_r8_cube.F \
 exch1_rl.F \
 exch1_rl_ad.F \
 exch1_rl_cube.F \
 exch1_rl_cube_ad.F \
 exch1_rs.F \
 exch1_rs_ad.F \
 exch1_rs_cube.F \
 exch1_rs_cube_ad.F \
 exch1_uv_r4_cube.F \
 exch1_uv_r8_cube.F \
 exch1_uv_rl_cube.F \
 exch1_uv_rs_cube.F \
 exch1_z_r4_cube.F \
 exch1_z_r8_cube.F \
 exch1_z_rl_cube.F \
 exch1_z_rs_cube.F \
 exch_3d_r4.F \
 exch_3d_r8.F \
 exch_cycle_ebl.F \
 exch_init.F \
 exch_r4_recv_get_x.F \
 exch_r4_recv_get_y.F \
 exch_r4_send_put_x.F \
 exch_r4_send_put_y.F \
 exch_r8_recv_get_x.F \
 exch_r8_recv_get_y.F \
 exch_r8_send_put_x.F \
 exch_r8_send_put_y.F \
 exch_rl_recv_get_x.F \
 exch_rl_recv_get_y.F \
 exch_rl_send_put_x.F \
 exch_rl_send_put_y.F \
 exch_rs_recv_get_x.F \
 exch_rs_recv_get_y.F \
 exch_rs_send_put_x.F \
 exch_rs_send_put_y.F \
 exch_s3d_r4.F \
 exch_s3d_r8.F \
 exch_s3d_rl.F \
 exch_s3d_rs.F \
 exch_sm_3d_r4.F \
 exch_sm_3d_r8.F \
 exch_sm_3d_rs.F \
 exch_tap_b.F \
 exch_tap_d.F \
 exch_uv_3d_r4.F \
 exch_uv_3d_r8.F \
 exch_uv_agrid_3d_r4.F \
 exch_uv_agrid_3d_r8.F \
 exch_uv_bgrid_3d_r4.F \
 exch_uv_bgrid_3d_r8.F \
 exch_uv_bgrid_3d_rl.F \
 exch_uv_bgrid_3d_rs.F \
 exch_uv_dgrid_3d_r4.F \
 exch_uv_dgrid_3d_r8.F \
 exch_uv_xy_r4.F \
 exch_uv_xy_r8.F \
 exch_uv_xyz_r4.F \
 exch_uv_xyz_r8.F \
 exch_xy_r4.F \
 exch_xy_r8.F \
 exch_xyz_r4.F \
 exch_xyz_r8.F \
 exch_z_3d_r4.F \
 exch_z_3d_r8.F \
 exch_z_3d_rl.F \
 exch_z_3d_rs.F \
 fool_the_compiler.F \
 gather_2d_r4.F \
 gather_2d_r8.F \
 gather_2d_wh_r4.F \
 gather_2d_wh_r8.F \
 gather_vec_r4.F \
 gather_vec_r8.F \
 gather_xz.F \
 gather_yz.F \
 global_max.F \
 global_sum.F \
 global_sum_singlecpu.F \
 global_sum_tile.F \
 global_sum_vector.F \
 ini_communication_patterns.F \
 ini_procs.F \
 ini_threading_environment.F \
 main.F \
 master_cpu_io.F \
 master_cpu_thread.F \
 mds_byteswapi4.F \
 mds_byteswapr4.F \
 mds_byteswapr8.F \
 mds_flush.F \
 mds_reclen.F \
 mdsfindunit.F \
 memsync.F \
 nml_change_syntax.F \
 nml_set_terminator.F \
 open_copy_data_file.F \
 print.F \
 reset_halo.F \
 scatter_2d_r4.F \
 scatter_2d_r8.F \
 scatter_2d_wh_r4.F \
 scatter_2d_wh_r8.F \
 scatter_vec_r4.F \
 scatter_vec_r8.F \
 scatter_xz.F \
 scatter_yz.F \
 stop_if_error.F \
 timers.F \
 utils.F \
 write_utils.F \
 add_walls2masks.F \
 calc_eddy_stress.F \
 calc_grid_angles.F \
 calc_gw.F \
 calc_oce_mxlayer.F \
 cg2d.F \
 cg2d_ex0.F \
 cg2d_sr.F \
 cg3d.F \
 cg3d_ex0.F \
 check_pickup.F \
 config_check.F \
 config_summary.F \
 diags_oceanic_surf_flux.F \
 do_statevars_diags.F \
 do_the_model_io.F \
 do_write_pickup.F \
 external_forcing.F \
 find_hyd_press_1d.F \
 gsw_teos10.F \
 ini_cartesian_grid.F \
 ini_cg3d.F \
 ini_cori.F \
 ini_curvilinear_grid.F \
 ini_cylinder_grid.F \
 ini_eos.F \
 ini_global_domain.F \
 ini_grid.F \
 ini_linear_phisurf.F \
 ini_local_grid.F \
 ini_masks_etc.F \
 ini_mnc_vars.F \
 ini_model_io.F \
 ini_nh_vars.F \
 ini_parms.F \
 ini_sigma_hfac.F \
 ini_spherical_polar_grid.F \
 ini_vertical_grid.F \
 initialise_fixed.F \
 load_ref_files.F \
 packages_boot.F \
 packages_check.F \
 packages_error_msg.F \
 packages_init_fixed.F \
 packages_print_msg.F \
 packages_readparms.F \
 packages_unused_msg.F \
 packages_write_pickup.F \
 plot_field.F \
 port_rand.F \
 post_cg3d.F \
 pre_cg3d.F \
 rotate_spherical_polar_grid.F \
 set_defaults.F \
 set_grid_factors.F \
 set_parms.F \
 set_ref_state.F \
 solve_uv_tridiago.F \
 the_model_main.F \
 timestep_wvel.F \
 tracers_iigw_correction.F \
 turnoff_model_io.F \
 write_grid.F \
 write_pickup.F \
 write_state.F
NON_AD_F90_SRC_FILES = 

# ... AD ...
adall: ad_exe_target
adtaf: ad_taf_output.f

ad_exe_target:
	@echo Update AD_CONFIG.h and make $(EXE_ADJ)
	@$(TOOLSDIR)/convert_cpp_cmd2defines "Adjoint version" -bAD_CONFIG_H -DALLOW_ADJOINT_RUN -UALLOW_TANGENTLINEAR_RUN > ad_config.template0
	@cmp ad_config.template0 AD_CONFIG.h || cat ad_config.template0 > AD_CONFIG.h
	@-rm -f ad_config.template0
	$(MAKE) -f $(MAKEFILE) $(EXE_ADJ)

ad_input_code.f: $(AD_FILES) $(AD_FLOW_FILES)
	@$(TOOLSDIR)/convert_cpp_cmd2defines "Adjoint version" -bAD_CONFIG_H -DALLOW_ADJOINT_RUN -UALLOW_TANGENTLINEAR_RUN > ad_config.template1
	cmp ad_config.template1 AD_CONFIG.h || cat ad_config.template1 > AD_CONFIG.h
	@-rm -f ad_config.template1
	@$(MAKE) -f $(MAKEFILE) $(F77_PP_SRC_FILES)
	cat $(AD_FLOW_FILES) $(AD_FILES) | sed -f $(TOOLSDIR)/remove_comments_sed > ad_input_code.f

# ... send 1 file to TAF ...
ad_taf_output.f: ad_input_code.f
	@-rm -f ad_input_code_ad.f ; echo ''
	$(TAF) $(AD_TAF_FLAGS) $(TAF_EXTRA) -fixed $(TAF_F77_FLAGS) ad_input_code.f
	ls -l ad_input_code_ad.f
	cat ad_input_code_ad.f | sed -f $(TOOLSDIR)/adjoint_sed > ad_taf_output.f

$(EXE_ADJ): $(SPECIAL_FILES) $(H_SRC_FILES) ad_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o) $(EMBEDDED_FILES)
	$(LINK) -o ${EXE_ADJ} $(FFLAGS) $(FOPTIM) ad_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o) $(LIBS)

adobj: ad_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o)

# ... FTL ...
ftlall: ftl_exe_target
ftltaf: ftl_taf_output.f

ftl_exe_target:
	@echo Update AD_CONFIG.h and make $(EXE_TLM)
	@$(TOOLSDIR)/convert_cpp_cmd2defines "TangLin version" -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -DALLOW_TANGENTLINEAR_RUN > ad_config.template0
	@cmp ad_config.template0 AD_CONFIG.h || cat ad_config.template0 > AD_CONFIG.h
	@-rm -f ad_config.template0
	$(MAKE) -f $(MAKEFILE) $(EXE_TLM)

ftl_input_code.f: $(AD_FILES) $(AD_FLOW_FILES)
	@$(TOOLSDIR)/convert_cpp_cmd2defines "TangLin version" -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -DALLOW_TANGENTLINEAR_RUN > ad_config.template1
	cmp ad_config.template1 AD_CONFIG.h || cat ad_config.template1 > AD_CONFIG.h
	@-rm -f ad_config.template1
	@$(MAKE) -f $(MAKEFILE) $(F77_PP_SRC_FILES)
	cat $(AD_FLOW_FILES) $(AD_FILES) | sed -f $(TOOLSDIR)/remove_comments_sed > ftl_input_code.f

# ... send 1 file to TAF ...
ftl_taf_output.f: ftl_input_code.f
	@-rm -f ftl_input_code_tl.f ; echo ''
	$(TAF) $(FTL_TAF_FLAGS) $(TAF_EXTRA) -fixed $(TAF_F77_FLAGS) ftl_input_code.f
	ls -l ftl_input_code_tl.f
	cat ftl_input_code_tl.f | sed -f $(TOOLSDIR)/adjoint_sed > ftl_taf_output.f

$(EXE_TLM): $(SPECIAL_FILES) $(H_SRC_FILES) ftl_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o) $(EMBEDDED_FILES)
	$(LINK) -o ${EXE_TLM} $(FFLAGS) $(FOPTIM) ftl_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o) $(LIBS)

ftlobj: ftl_taf_output.o $(NON_AD_F90_SRC_FILES:.F90=.o) $(NON_AD_F77_SRC_FILES:.F=.o) $(C_SRC_FILES:.c=.o)

# ... SVD ...
svdtaf: ad_taf_output.f ftl_taf_output.f
	@echo "--->>> Only ran TAF to generate SVD code!    <<<---"
	@echo "--->>> Do make svdall afterwards to compile. <<<---"
svdall: svd_touch svd_taf

svd_taf: $(OBJFILES)
	$(LINK) -o mitgcmuv_svd $(FFLAGS) $(FOPTIM) $(OBJFILES) ad_taf_output.o ftl_taf_output.o $(LIBS)

	@echo "--->>> Only COMPILE svd code! <<<---"
	@echo "--->>> Assumes you previously <<<---"
	@echo "--->>> did make svdtaf        <<<---"

svd_touch:
	@echo "--->>> Only COMPILE svd code! <<<---"
	@echo "--->>> Assumes you previously <<<---"
	@echo "--->>> did make svdtaf        <<<---"
	touch ad_taf_output.f ftl_taf_output.f
	$(FC) $(FFLAGS) $(FOPTIM) -c ad_taf_output.f
	$(FC) $(FFLAGS) $(FOPTIM) -c ftl_taf_output.f
	@$(TOOLSDIR)/convert_cpp_cmd2defines "SVD version" -bAD_CONFIG_H -UALLOW_ADJOINT_RUN -DALLOW_TANGENTLINEAR_RUN > ad_config.template
	cmp ad_config.template AD_CONFIG.h || cat ad_config.template > AD_CONFIG.h
	@-rm -f ad_config.template

ini_masks_etc.o: ini_masks_etc.f
	$(FC) $(FFLAGS) $(NOOPTFLAGS) -c $<
# This section creates symbolic links


#  These files are linked from /home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code ( aka /home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code )
 DIAGNOSTICS_SIZE.h SIZE.h :
	if [ ! -L $@ ]; then $(LN) /home/hmh85/scratch/MITgcm/verification/vortexSphere_Williamson_TC2/code/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/debug ( aka /home/hmh85/scratch/MITgcm/pkg/debug )
 DEBUG_OPTIONS.h chksum_tiled.F debug_call.F debug_cs_corner_uv.F debug_enter.F debug_fld_stats_rl.F debug_fld_stats_rs.F debug_leave.F debug_msg.F debug_stats_rl.F debug_stats_rs.F fill_in_corners_rl.F write_fullarray_rl.F write_fullarray_rs.F debug_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/debug/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/diagnostics ( aka /home/hmh85/scratch/MITgcm/pkg/diagnostics )
 DIAGNOSTICS.h DIAGNOSTICS_CALC.h DIAGNOSTICS_P2SHARE.h DIAGSTATS_REGIONS.h DIAG_OPTIONS.h diag_calc_psivel.F diag_cg2d.F diag_vegtile_fill.F diagnostics_addtolist.F diagnostics_calc_phivel.F diagnostics_check.F diagnostics_clear.F diagnostics_fill.F diagnostics_fill_field.F diagnostics_fill_rs.F diagnostics_fill_state.F diagnostics_fract_fill.F diagnostics_ini_io.F diagnostics_init_early.F diagnostics_init_fixed.F diagnostics_init_varia.F diagnostics_interp_p2p.F diagnostics_interp_vert.F diagnostics_is_on.F diagnostics_list_check.F diagnostics_main_init.F diagnostics_mnc_out.F diagnostics_out.F diagnostics_read_pickup.F diagnostics_readparms.F diagnostics_scale_fill.F diagnostics_scale_fill_rs.F diagnostics_set_calc.F diagnostics_set_levels.F diagnostics_set_pointers.F diagnostics_setdiag.F diagnostics_status_error.F diagnostics_sum_levels.F diagnostics_summary.F diagnostics_switch_onoff.F diagnostics_utils.F diagnostics_write.F diagnostics_write_adj.F diagnostics_write_pickup.F diagstats_ascii_out.F diagstats_calc.F diagstats_clear.F diagstats_close_io.F diagstats_fill.F diagstats_global.F diagstats_ini_io.F diagstats_local.F diagstats_mnc_out.F diagstats_others_calc.F diagstats_output.F diagstats_set_pointers.F diagstats_set_regions.F diagstats_setdiag.F diagnostics_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/diagnostics/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/generic_advdiff ( aka /home/hmh85/scratch/MITgcm/pkg/generic_advdiff )
 GAD.h GAD_FLUX_LIMITER.h GAD_OPTIONS.h GAD_SOM_VARS.h gad_ad_check_lev1_dir.h gad_ad_check_lev2_dir.h gad_ad_check_lev3_dir.h gad_ad_check_lev4_dir.h gad_advection.F gad_advscheme.F gad_biharm_r.F gad_biharm_x.F gad_biharm_y.F gad_c2_adv_r.F gad_c2_adv_x.F gad_c2_adv_y.F gad_c2_impl_r.F gad_c4_adv_r.F gad_c4_adv_x.F gad_c4_adv_y.F gad_calc_rhs.F gad_check.F gad_del2.F gad_diagnostics_init.F gad_diagnostics_state.F gad_diff_r.F gad_diff_x.F gad_diff_y.F gad_dst2u1_adv_r.F gad_dst2u1_adv_x.F gad_dst2u1_adv_y.F gad_dst2u1_impl_r.F gad_dst3_adv_r.F gad_dst3_adv_x.F gad_dst3_adv_y.F gad_dst3fl_adv_r.F gad_dst3fl_adv_x.F gad_dst3fl_adv_y.F gad_dst3fl_impl_r.F gad_exch_som.F gad_fluxlimit_adv_r.F gad_fluxlimit_adv_x.F gad_fluxlimit_adv_y.F gad_fluxlimit_impl_r.F gad_grad_x.F gad_grad_y.F gad_implicit_r.F gad_init_fixed.F gad_init_varia.F gad_os7mp_adv_r.F gad_os7mp_adv_x.F gad_os7mp_adv_y.F gad_osc_hat_r.F gad_osc_hat_x.F gad_osc_hat_y.F gad_osc_mul_r.F gad_osc_mul_x.F gad_osc_mul_y.F gad_plm_fun.F gad_ppm_adv_r.F gad_ppm_adv_x.F gad_ppm_adv_y.F gad_ppm_flx_r.F gad_ppm_flx_x.F gad_ppm_flx_y.F gad_ppm_fun.F gad_ppm_hat_r.F gad_ppm_hat_x.F gad_ppm_hat_y.F gad_ppm_p3e_r.F gad_ppm_p3e_x.F gad_ppm_p3e_y.F gad_pqm_adv_r.F gad_pqm_adv_x.F gad_pqm_adv_y.F gad_pqm_flx_r.F gad_pqm_flx_x.F gad_pqm_flx_y.F gad_pqm_fun.F gad_pqm_hat_r.F gad_pqm_hat_x.F gad_pqm_hat_y.F gad_pqm_p5e_r.F gad_pqm_p5e_x.F gad_pqm_p5e_y.F gad_read_pickup.F gad_som_adv_r.F gad_som_adv_x.F gad_som_adv_y.F gad_som_advect.F gad_som_exchanges.F gad_som_fill_cs_corner.F gad_som_lim_r.F gad_som_prep_cs_corner.F gad_u3_adv_r.F gad_u3_adv_x.F gad_u3_adv_y.F gad_u3c4_impl_r.F gad_write_pickup.F salt_fill.F :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/generic_advdiff/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/mdsio ( aka /home/hmh85/scratch/MITgcm/pkg/mdsio )
 MDSIO_BUFF_3D.h MDSIO_BUFF_WH.h MDSIO_OPTIONS.h mdsio_buffertorl.F mdsio_buffertors.F mdsio_check4file.F mdsio_facef_read.F mdsio_gl.F mdsio_gl_slice.F mdsio_pass_r4torl.F mdsio_pass_r4tors.F mdsio_pass_r8torl.F mdsio_pass_r8tors.F mdsio_rd_rec_rl.F mdsio_rd_rec_rs.F mdsio_read_field.F mdsio_read_meta.F mdsio_read_section.F mdsio_read_tape.F mdsio_read_whalos.F mdsio_readvec_loc.F mdsio_rw_field.F mdsio_rw_slice.F mdsio_seg4torl.F mdsio_seg4tors.F mdsio_seg8torl.F mdsio_seg8tors.F mdsio_segxtorx_2d.F mdsio_wr_metafiles.F mdsio_wr_rec_rl.F mdsio_wr_rec_rs.F mdsio_write_field.F mdsio_write_meta.F mdsio_write_section.F mdsio_write_tape.F mdsio_write_whalos.F mdsio_writelocal.F mdsio_writevec_loc.F mdsio_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/mdsio/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/mom_common ( aka /home/hmh85/scratch/MITgcm/pkg/mom_common )
 MOM_COMMON_OPTIONS.h MOM_VISC.h mom_calc_3d_strain.F mom_calc_absvort3.F mom_calc_hdiv.F mom_calc_hfacz.F mom_calc_ke.F mom_calc_relvort3.F mom_calc_smag_3d.F mom_calc_strain.F mom_calc_tension.F mom_calc_visc.F mom_diagnostics_init.F mom_hdissip.F mom_init_fixed.F mom_quasihydrostatic.F mom_u_botdrag_coeff.F mom_u_coriolis_nh.F mom_u_implicit_r.F mom_u_metric_nh.F mom_u_rviscflux.F mom_u_sidedrag.F mom_uv_smag_3d.F mom_v_botdrag_coeff.F mom_v_coriolis_nh.F mom_v_implicit_r.F mom_v_metric_nh.F mom_v_rviscflux.F mom_v_sidedrag.F mom_visc_qgl_limit.F mom_visc_qgl_stretch.F mom_w_coriolis_nh.F mom_w_metric_nh.F mom_w_sidedrag.F mom_w_smag_3d.F :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/mom_common/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/mom_fluxform ( aka /home/hmh85/scratch/MITgcm/pkg/mom_fluxform )
 MOM_FLUXFORM.h MOM_FLUXFORM_OPTIONS.h mom_calc_rtrans.F mom_fluxform.F mom_u_adv_uu.F mom_u_adv_vu.F mom_u_adv_wu.F mom_u_coriolis.F mom_u_del2u.F mom_u_metric_cylinder.F mom_u_metric_sphere.F mom_u_xviscflux.F mom_u_yviscflux.F mom_uv_boundary.F mom_v_adv_uv.F mom_v_adv_vv.F mom_v_adv_wv.F mom_v_coriolis.F mom_v_del2v.F mom_v_metric_cylinder.F mom_v_metric_sphere.F mom_v_xviscflux.F mom_v_yviscflux.F :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/mom_fluxform/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/mom_vecinv ( aka /home/hmh85/scratch/MITgcm/pkg/mom_vecinv )
 MOM_VECINV_OPTIONS.h mom_vecinv.F mom_vi_coriolis.F mom_vi_del2uv.F mom_vi_hdissip.F mom_vi_u_coriolis.F mom_vi_u_coriolis_c4.F mom_vi_u_grad_ke.F mom_vi_u_vertshear.F mom_vi_v_coriolis.F mom_vi_v_coriolis_c4.F mom_vi_v_grad_ke.F mom_vi_v_vertshear.F :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/mom_vecinv/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/monitor ( aka /home/hmh85/scratch/MITgcm/pkg/monitor )
 MONITOR.h MONITOR_OPTIONS.h mon_advcfl.F mon_advcflw.F mon_advcflw2.F mon_calc_advcfl.F mon_calc_stats_rl.F mon_calc_stats_rs.F mon_init.F mon_ke.F mon_out.F mon_printstats_rl.F mon_printstats_rs.F mon_set_iounit.F mon_set_pref.F mon_solution.F mon_stats_latbnd_rl.F mon_stats_rl.F mon_stats_rs.F mon_surfcor.F mon_vort3.F mon_writestats_rl.F mon_writestats_rs.F monitor.F monitor_ad.F monitor_g.F monitor_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/monitor/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/mypackage ( aka /home/hmh85/scratch/MITgcm/pkg/mypackage )
 MYPACKAGE.h MYPACKAGE_OPTIONS.h mypackage_calc_rhs.F mypackage_check.F mypackage_diagnostics_init.F mypackage_diagnostics_state.F mypackage_init_fixed.F mypackage_init_varia.F mypackage_mnc_init.F mypackage_read_pickup.F mypackage_readparms.F mypackage_tendency_apply.F mypackage_write_pickup.F mypackage_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/mypackage/$@ $@; fi

#  These files are linked from $(ROOTDIR)/pkg/rw ( aka /home/hmh85/scratch/MITgcm/pkg/rw )
 RW_MFLDS.h RW_OPTIONS.h get_write_global_fld.F read_fld_xy_rl.F read_fld_xy_rs.F read_fld_xyz_rl.F read_fld_xyz_rs.F read_glvec_rl.F read_glvec_rs.F read_mflds.F read_rec.F rw_get_suffix.F set_write_global_fld.F write_fld_3d_rl.F write_fld_3d_rs.F write_fld_s3d_rl.F write_fld_s3d_rs.F write_fld_xy_rl.F write_fld_xy_rs.F write_fld_xyz_rl.F write_fld_xyz_rs.F write_glvec_rl.F write_glvec_rs.F write_local_rl.F write_local_rs.F write_rec.F rw_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/pkg/rw/$@ $@; fi

#  These files are linked from $(ROOTDIR)/eesupp/src ( aka /home/hmh85/scratch/MITgcm/eesupp/src )
 all_proc_die.F bar2.F bar_check.F barrier.F check_threads.F comm_stats.F cumulsum_z_tile.F diff_phase_multiple.F different_multiple.F eeboot.F eeboot_minimal.F eedata_example.F eedie.F eeintro_msg.F eeset_parms.F eewrite_eeenv.F exch0_r4.F exch0_r8.F exch0_rl.F exch0_rs.F exch1_bg_r4_cube.F exch1_bg_r8_cube.F exch1_bg_rl_cube.F exch1_bg_rs_cube.F exch1_r4.F exch1_r4_cube.F exch1_r8.F exch1_r8_cube.F exch1_rl.F exch1_rl_ad.F exch1_rl_cube.F exch1_rl_cube_ad.F exch1_rs.F exch1_rs_ad.F exch1_rs_cube.F exch1_rs_cube_ad.F exch1_uv_r4_cube.F exch1_uv_r8_cube.F exch1_uv_rl_cube.F exch1_uv_rs_cube.F exch1_z_r4_cube.F exch1_z_r8_cube.F exch1_z_rl_cube.F exch1_z_rs_cube.F exch_3d_r4.F exch_3d_r8.F exch_3d_rl.F exch_3d_rs.F exch_cycle_ebl.F exch_init.F exch_r4_recv_get_x.F exch_r4_recv_get_y.F exch_r4_send_put_x.F exch_r4_send_put_y.F exch_r8_recv_get_x.F exch_r8_recv_get_y.F exch_r8_send_put_x.F exch_r8_send_put_y.F exch_rl_recv_get_x.F exch_rl_recv_get_y.F exch_rl_send_put_x.F exch_rl_send_put_y.F exch_rs_recv_get_x.F exch_rs_recv_get_y.F exch_rs_send_put_x.F exch_rs_send_put_y.F exch_s3d_r4.F exch_s3d_r8.F exch_s3d_rl.F exch_s3d_rs.F exch_sm_3d_r4.F exch_sm_3d_r8.F exch_sm_3d_rl.F exch_sm_3d_rs.F exch_tap_b.F exch_tap_d.F exch_uv_3d_r4.F exch_uv_3d_r8.F exch_uv_3d_rl.F exch_uv_3d_rs.F exch_uv_agrid_3d_r4.F exch_uv_agrid_3d_r8.F exch_uv_agrid_3d_rl.F exch_uv_agrid_3d_rs.F exch_uv_bgrid_3d_r4.F exch_uv_bgrid_3d_r8.F exch_uv_bgrid_3d_rl.F exch_uv_bgrid_3d_rs.F exch_uv_dgrid_3d_r4.F exch_uv_dgrid_3d_r8.F exch_uv_dgrid_3d_rl.F exch_uv_dgrid_3d_rs.F exch_uv_xy_r4.F exch_uv_xy_r8.F exch_uv_xy_rl.F exch_uv_xy_rs.F exch_uv_xyz_r4.F exch_uv_xyz_r8.F exch_uv_xyz_rl.F exch_uv_xyz_rs.F exch_xy_r4.F exch_xy_r8.F exch_xy_rl.F exch_xy_rs.F exch_xyz_r4.F exch_xyz_r8.F exch_xyz_rl.F exch_xyz_rs.F exch_z_3d_r4.F exch_z_3d_r8.F exch_z_3d_rl.F exch_z_3d_rs.F fill_cs_corner_ag_rl.F fill_cs_corner_tr_rl.F fill_cs_corner_uv_rl.F fill_cs_corner_uv_rs.F fool_the_compiler.F gather_2d_r4.F gather_2d_r8.F gather_2d_wh_r4.F gather_2d_wh_r8.F gather_vec_r4.F gather_vec_r8.F gather_xz.F gather_yz.F get_periodic_interval.F global_max.F global_sum.F global_sum_singlecpu.F global_sum_tile.F global_sum_vector.F gsl_ieee_env.c ini_communication_patterns.F ini_procs.F ini_threading_environment.F lef_zero.F main.F master_cpu_io.F master_cpu_thread.F mds_byteswapi4.F mds_byteswapr4.F mds_byteswapr8.F mds_flush.F mds_reclen.F mdsfindunit.F memsync.F nml_change_syntax.F nml_set_terminator.F open_copy_data_file.F print.F ptwrapper.c reset_halo.F scatter_2d_r4.F scatter_2d_r8.F scatter_2d_wh_r4.F scatter_2d_wh_r8.F scatter_vec_r4.F scatter_vec_r8.F scatter_xz.F scatter_yz.F setdir.c setrlstk.c sigreg.c stop_if_error.F tim.c timer_stats.c timers.F utils.F write_utils.F eesupp.flow exch1_ad.flow :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/eesupp/src/$@ $@; fi

#  These files are linked from $(ROOTDIR)/model/src ( aka /home/hmh85/scratch/MITgcm/model/src )
 adams_bashforth2.F adams_bashforth3.F add_walls2masks.F apply_forcing.F calc_3d_diffusivity.F calc_adv_flow.F calc_div_ghat.F calc_eddy_stress.F calc_grad_phi_fv.F calc_grad_phi_hyd.F calc_grad_phi_surf.F calc_grid_angles.F calc_gw.F calc_ivdc.F calc_oce_mxlayer.F calc_phi_hyd.F calc_r_star.F calc_surf_dr.F calc_viscosity.F calc_wsurf_tr.F cg2d.F cg2d_ex0.F cg2d_nsa.F cg2d_sr.F cg3d.F cg3d_ex0.F check_pickup.F config_check.F config_summary.F convective_adjustment.F convective_adjustment_ini.F convective_weights.F convectively_mixtracer.F correction_step.F cycle_ab_tracer.F cycle_tracer.F diags_oceanic_surf_flux.F diags_phi_hyd.F diags_phi_rlow.F diags_rho.F diags_sound_speed.F do_atmospheric_phys.F do_fields_blocking_exchanges.F do_oceanic_phys.F do_stagger_fields_exchanges.F do_statevars_diags.F do_the_model_io.F do_write_pickup.F dynamics.F external_fields_load.F external_forcing.F external_forcing_surf.F find_alpha.F find_hyd_press_1d.F find_rho.F forcing_surf_relax.F forward_step.F freesurf_rescale_g.F freeze_surface.F grad_sigma.F gsw_teos10.F impldiff.F ini_cartesian_grid.F ini_cg2d.F ini_cg3d.F ini_cori.F ini_curvilinear_grid.F ini_cylinder_grid.F ini_depths.F ini_dynvars.F ini_eos.F ini_ffields.F ini_fields.F ini_forcing.F ini_global_domain.F ini_grid.F ini_linear_phisurf.F ini_local_grid.F ini_masks_etc.F ini_mixing.F ini_mnc_vars.F ini_model_io.F ini_nh_fields.F ini_nh_vars.F ini_nlfs_vars.F ini_p_ground.F ini_parms.F ini_pressure.F ini_psurf.F ini_salt.F ini_sigma_hfac.F ini_spherical_polar_grid.F ini_theta.F ini_vel.F ini_vertical_grid.F initialise_fixed.F initialise_varia.F integr_continuity.F integrate_for_w.F load_fields_driver.F load_grid_spacing.F load_ref_files.F main_do_loop.F momentum_correction_step.F packages_boot.F packages_check.F packages_error_msg.F packages_init_fixed.F packages_init_variables.F packages_print_msg.F packages_readparms.F packages_unused_msg.F packages_write_pickup.F plot_field.F port_rand.F post_cg3d.F pre_cg3d.F pressure_for_eos.F read_pickup.F remove_mean.F reset_nlfs_vars.F rotate_spherical_polar_grid.F rotate_uv2en.F salt_integrate.F seawater.F set_defaults.F set_grid_factors.F set_parms.F set_ref_state.F solve_for_pressure.F solve_pentadiagonal.F solve_tridiagonal.F solve_uv_tridiago.F state_summary.F swfrac.F taueddy_init_varia.F taueddy_tendency_apply.F temp_integrate.F the_main_loop.F the_model_main.F thermodynamics.F timestep.F timestep_tracer.F timestep_wvel.F tracers_correction_step.F tracers_iigw_correction.F turnoff_model_io.F update_cg2d.F update_etah.F update_etaws.F update_masks_etc.F update_r_star.F update_sigma.F update_surf_dr.F write_grid.F write_pickup.F write_state.F :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/model/src/$@ $@; fi

#  These files are linked from /apps/sw/mpi/openmpi/4.1.4-gcc-9.1.0-slurm-18.08.6/include ( aka /apps/sw/mpi/openmpi/4.1.4-gcc-9.1.0-slurm-18.08.6/include )
 mpi-ext.h mpi.h mpi_portable_platform.h mpif-c-constants-decl.h mpif-config.h mpif-constants.h mpif-ext.h mpif-externals.h mpif-handles.h mpif-io-constants.h mpif-io-handles.h mpif-sentinels.h mpif-sizeof.h mpif.h pshmem.h pshmemx.h shmem-compat.h shmem.h shmemx.h :
	if [ ! -L $@ ]; then $(LN) /apps/sw/mpi/openmpi/4.1.4-gcc-9.1.0-slurm-18.08.6/include/$@ $@; fi

#  These files are linked from $(ROOTDIR)/eesupp/inc ( aka /home/hmh85/scratch/MITgcm/eesupp/inc )
 BAR2.h BARRIER.h CPP_EEMACROS.h CPP_EEOPTIONS.h CUMULSUM.h DEF_IN_MAKEFILE.h EEBUFF_SCPU.h EEPARAMS.h EESUPPORT.h EXCH.h GLOBAL_MAX.h GLOBAL_SUM.h MAIN_PDIRECTIVES1.h MAIN_PDIRECTIVES2.h MPI_INFO.h SIGREG.h :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/eesupp/inc/$@ $@; fi

#  These files are linked from $(ROOTDIR)/model/inc ( aka /home/hmh85/scratch/MITgcm/model/inc )
 CG2D.h CG3D.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h EOS.h FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h NH_VARS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h RESTART.h SET_GRID.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h :
	if [ ! -L $@ ]; then $(LN) $(ROOTDIR)/model/inc/$@ $@; fi


# DO NOT DELETE


adams_bashforth2.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
adams_bashforth2.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
adams_bashforth3.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
adams_bashforth3.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
add_walls2masks.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
add_walls2masks.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
add_walls2masks.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
add_walls2masks.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
add_walls2masks.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
add_walls2masks.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
add_walls2masks.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
add_walls2masks.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
add_walls2masks.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
add_walls2masks.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
add_walls2masks.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
add_walls2masks.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
all_proc_die.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
all_proc_die.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
apply_forcing.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
apply_forcing.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
apply_forcing.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
apply_forcing.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
apply_forcing.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
apply_forcing.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
apply_forcing.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
apply_forcing.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
apply_forcing.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
apply_forcing.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
apply_forcing.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
apply_forcing.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
bar2.f: /usr/include/stdc-predef.h BAR2.h CPP_EEMACROS.h CPP_EEOPTIONS.h
bar2.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
bar_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
bar_check.f: EEPARAMS.h SIZE.h
barrier.f: /usr/include/stdc-predef.h BARRIER.h CPP_EEMACROS.h CPP_EEOPTIONS.h
barrier.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
calc_3d_diffusivity.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
calc_3d_diffusivity.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
calc_3d_diffusivity.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
calc_3d_diffusivity.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
calc_3d_diffusivity.f: FCORI_MACROS.h GAD.h GRID.h GRID_MACROS.h
calc_3d_diffusivity.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_3d_diffusivity.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
calc_3d_diffusivity.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
calc_3d_diffusivity.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
calc_3d_diffusivity.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_3d_diffusivity.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
calc_3d_diffusivity.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
calc_3d_diffusivity.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
calc_3d_diffusivity.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
calc_3d_diffusivity.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
calc_adv_flow.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_adv_flow.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_adv_flow.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_adv_flow.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_adv_flow.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_adv_flow.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
calc_adv_flow.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_adv_flow.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_adv_flow.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_adv_flow.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
calc_adv_flow.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
calc_adv_flow.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
calc_div_ghat.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_div_ghat.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_div_ghat.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
calc_div_ghat.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_div_ghat.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_div_ghat.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
calc_div_ghat.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
calc_div_ghat.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
calc_div_ghat.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_div_ghat.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_div_ghat.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
calc_eddy_stress.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_eddy_stress.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_eddy_stress.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_eddy_stress.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
calc_eddy_stress.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_eddy_stress.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
calc_eddy_stress.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_eddy_stress.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_eddy_stress.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_eddy_stress.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_eddy_stress.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_eddy_stress.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
calc_eddy_stress.f: YC_MACROS.h
calc_grad_phi_fv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_grad_phi_fv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_grad_phi_fv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_grad_phi_fv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
calc_grad_phi_fv.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_grad_phi_fv.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
calc_grad_phi_fv.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_grad_phi_fv.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_grad_phi_fv.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_grad_phi_fv.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_grad_phi_fv.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
calc_grad_phi_fv.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
calc_grad_phi_fv.f: YC_MACROS.h
calc_grad_phi_hyd.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_grad_phi_hyd.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_grad_phi_hyd.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_grad_phi_hyd.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
calc_grad_phi_hyd.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
calc_grad_phi_hyd.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
calc_grad_phi_hyd.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_grad_phi_hyd.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_grad_phi_hyd.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_grad_phi_hyd.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_grad_phi_hyd.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_grad_phi_hyd.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_grad_phi_hyd.f: XC_MACROS.h YC_MACROS.h
calc_grad_phi_surf.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
calc_grad_phi_surf.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
calc_grad_phi_surf.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
calc_grad_phi_surf.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
calc_grad_phi_surf.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
calc_grad_phi_surf.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
calc_grad_phi_surf.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_grad_phi_surf.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_grad_phi_surf.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_grad_phi_surf.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_grad_phi_surf.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_grad_phi_surf.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_grad_phi_surf.f: XC_MACROS.h YC_MACROS.h
calc_grid_angles.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_grid_angles.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_grid_angles.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_grid_angles.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
calc_grid_angles.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_grid_angles.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
calc_grid_angles.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_grid_angles.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_grid_angles.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_grid_angles.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_grid_angles.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_grid_angles.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
calc_grid_angles.f: YC_MACROS.h
calc_gw.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_gw.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
calc_gw.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
calc_gw.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
calc_gw.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
calc_gw.f: MOM_COMMON_OPTIONS.h MOM_VISC.h NH_VARS.h PACKAGES_CONFIG.h
calc_gw.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_gw.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_gw.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_gw.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
calc_gw.f: RECIP_HFACW_MACROS.h RESTART.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
calc_gw.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
calc_ivdc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_ivdc.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
calc_ivdc.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
calc_ivdc.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
calc_ivdc.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
calc_ivdc.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_ivdc.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_ivdc.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_ivdc.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
calc_ivdc.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_ivdc.f: XC_MACROS.h YC_MACROS.h
calc_oce_mxlayer.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_oce_mxlayer.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_oce_mxlayer.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_oce_mxlayer.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
calc_oce_mxlayer.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_oce_mxlayer.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
calc_oce_mxlayer.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_oce_mxlayer.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_oce_mxlayer.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_oce_mxlayer.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_oce_mxlayer.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
calc_oce_mxlayer.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
calc_oce_mxlayer.f: YC_MACROS.h
calc_phi_hyd.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_phi_hyd.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_phi_hyd.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
calc_phi_hyd.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_phi_hyd.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_phi_hyd.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
calc_phi_hyd.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
calc_phi_hyd.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
calc_phi_hyd.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
calc_phi_hyd.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
calc_phi_hyd.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_phi_hyd.f: XC_MACROS.h YC_MACROS.h
calc_r_star.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_r_star.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_r_star.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_r_star.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_r_star.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_r_star.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
calc_r_star.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
calc_r_star.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
calc_r_star.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
calc_r_star.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
calc_r_star.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_r_star.f: XC_MACROS.h YC_MACROS.h
calc_surf_dr.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_surf_dr.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_surf_dr.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_surf_dr.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_surf_dr.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_surf_dr.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
calc_surf_dr.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
calc_surf_dr.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
calc_surf_dr.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
calc_surf_dr.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
calc_surf_dr.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
calc_surf_dr.f: XC_MACROS.h YC_MACROS.h
calc_viscosity.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_viscosity.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_viscosity.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_viscosity.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
calc_viscosity.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
calc_viscosity.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
calc_viscosity.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
calc_viscosity.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
calc_viscosity.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
calc_viscosity.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
calc_viscosity.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
calc_viscosity.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
calc_wsurf_tr.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
calc_wsurf_tr.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
calc_wsurf_tr.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
calc_wsurf_tr.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
calc_wsurf_tr.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
calc_wsurf_tr.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
calc_wsurf_tr.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
calc_wsurf_tr.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
calc_wsurf_tr.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
calc_wsurf_tr.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
calc_wsurf_tr.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
cg2d.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg2d.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
cg2d_ex0.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg2d_ex0.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
cg2d_nsa.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg2d_nsa.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
cg2d_sr.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg2d_sr.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
cg3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg3d.f: CPP_OPTIONS.h EEPARAMS.h SIZE.h
cg3d_ex0.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cg3d_ex0.f: CPP_OPTIONS.h EEPARAMS.h SIZE.h
check_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
check_pickup.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h RESTART.h SIZE.h
check_threads.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
check_threads.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
chksum_tiled.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
chksum_tiled.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
chksum_tiled.f: SIZE.h
comm_stats.f: /usr/include/stdc-predef.h BAR2.h CPP_EEMACROS.h CPP_EEOPTIONS.h
comm_stats.f: EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h mpif.h
config_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
config_check.f: CPP_OPTIONS.h EEPARAMS.h MOM_COMMON_OPTIONS.h MOM_VISC.h
config_check.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
config_summary.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
config_summary.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
config_summary.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
config_summary.f: DYU_MACROS.h EEPARAMS.h EOS.h FCORI_MACROS.h GRID.h
config_summary.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
config_summary.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
config_summary.f: MOM_VISC.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
config_summary.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
config_summary.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
config_summary.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
config_summary.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
config_summary.f: RECIP_HFACW_MACROS.h SET_GRID.h SIZE.h TANPHIATU_MACROS.h
config_summary.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
convective_adjustment.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
convective_adjustment.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
convective_adjustment.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
convective_adjustment.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
convective_adjustment.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
convective_adjustment.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
convective_adjustment.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
convective_adjustment.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
convective_adjustment.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
convective_adjustment.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
convective_adjustment.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
convective_adjustment.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
convective_adjustment.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
convective_adjustment.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
convective_adjustment.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
convective_adjustment_ini.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
convective_adjustment_ini.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
convective_adjustment_ini.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
convective_adjustment_ini.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
convective_adjustment_ini.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
convective_adjustment_ini.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
convective_adjustment_ini.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
convective_adjustment_ini.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
convective_adjustment_ini.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
convective_adjustment_ini.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
convective_adjustment_ini.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
convective_adjustment_ini.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
convective_adjustment_ini.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
convective_adjustment_ini.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
convective_adjustment_ini.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
convective_adjustment_ini.f: YC_MACROS.h
convective_weights.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
convective_weights.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
convective_weights.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
convective_weights.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
convective_weights.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
convective_weights.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
convective_weights.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
convective_weights.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
convective_weights.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
convective_weights.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
convective_weights.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
convective_weights.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
convective_weights.f: YC_MACROS.h
convectively_mixtracer.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
convectively_mixtracer.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h SIZE.h
correction_step.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
correction_step.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
correction_step.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
correction_step.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
correction_step.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
correction_step.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
correction_step.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
correction_step.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
correction_step.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
correction_step.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
correction_step.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
correction_step.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
correction_step.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
cumulsum_z_tile.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cumulsum_z_tile.f: CUMULSUM.h EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h
cumulsum_z_tile.f: mpif.h
cycle_ab_tracer.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cycle_ab_tracer.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
cycle_tracer.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
cycle_tracer.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
debug_call.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_call.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
debug_call.f: SIZE.h
debug_cs_corner_uv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
debug_cs_corner_uv.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h
debug_cs_corner_uv.f: PACKAGES_CONFIG.h SIZE.h
debug_enter.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_enter.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
debug_enter.f: SIZE.h
debug_fld_stats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
debug_fld_stats_rl.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h
debug_fld_stats_rl.f: PACKAGES_CONFIG.h SIZE.h
debug_fld_stats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
debug_fld_stats_rs.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h
debug_fld_stats_rs.f: PACKAGES_CONFIG.h SIZE.h
debug_leave.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_leave.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
debug_leave.f: SIZE.h
debug_msg.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_msg.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
debug_stats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_stats_rl.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
debug_stats_rl.f: SIZE.h
debug_stats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
debug_stats_rs.f: CPP_OPTIONS.h DEBUG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
debug_stats_rs.f: SIZE.h
diag_calc_psivel.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diag_calc_psivel.f: CPP_OPTIONS.h DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
diag_calc_psivel.f: PARAMS.h SIZE.h
diag_cg2d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diag_cg2d.f: CPP_OPTIONS.h DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
diag_cg2d.f: PARAMS.h SIZE.h
diag_vegtile_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diag_vegtile_fill.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diag_vegtile_fill.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagnostics_addtolist.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_addtolist.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_addtolist.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_addtolist.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_calc_phivel.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_calc_phivel.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_calc_phivel.f: DIAGNOSTICS_CALC.h DIAGNOSTICS_SIZE.h
diagnostics_calc_phivel.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
diagnostics_calc_phivel.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
diagnostics_calc_phivel.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
diagnostics_calc_phivel.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diagnostics_calc_phivel.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diagnostics_calc_phivel.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagnostics_calc_phivel.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_calc_phivel.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
diagnostics_calc_phivel.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diagnostics_calc_phivel.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_calc_phivel.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_calc_phivel.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_calc_phivel.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagnostics_calc_phivel.f: YC_MACROS.h
diagnostics_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_check.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_P2SHARE.h
diagnostics_check.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h DXC_MACROS.h
diagnostics_check.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diagnostics_check.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diagnostics_check.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
diagnostics_check.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_check.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagnostics_check.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_check.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagnostics_check.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_check.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_check.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_check.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagnostics_check.f: YC_MACROS.h
diagnostics_clear.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_clear.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagnostics_clear.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagnostics_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_fill.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagnostics_fill.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagnostics_fill_field.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_fill_field.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_fill_field.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h DXC_MACROS.h
diagnostics_fill_field.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diagnostics_fill_field.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diagnostics_fill_field.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
diagnostics_fill_field.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_fill_field.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
diagnostics_fill_field.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_fill_field.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
diagnostics_fill_field.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diagnostics_fill_field.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_fill_field.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_fill_field.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_fill_field.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagnostics_fill_field.f: YC_MACROS.h
diagnostics_fill_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_fill_rs.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_fill_rs.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_fill_rs.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_fill_state.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_fill_state.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAG_OPTIONS.h
diagnostics_fill_state.f: DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
diagnostics_fill_state.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
diagnostics_fill_state.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GRID.h
diagnostics_fill_state.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diagnostics_fill_state.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diagnostics_fill_state.f: NH_VARS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagnostics_fill_state.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_fill_state.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
diagnostics_fill_state.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diagnostics_fill_state.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_fill_state.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_fill_state.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_fill_state.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
diagnostics_fill_state.f: XC_MACROS.h YC_MACROS.h
diagnostics_fract_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_fract_fill.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_fract_fill.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_fract_fill.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_ini_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_ini_io.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_ini_io.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_ini_io.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_init_early.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_init_early.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_init_early.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_init_early.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_init_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_init_fixed.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_init_fixed.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_init_fixed.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_init_varia.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_init_varia.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_init_varia.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_init_varia.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_interp_p2p.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_interp_p2p.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAG_OPTIONS.h
diagnostics_interp_p2p.f: PACKAGES_CONFIG.h
diagnostics_interp_vert.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_interp_vert.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_interp_vert.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h DXC_MACROS.h
diagnostics_interp_vert.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diagnostics_interp_vert.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diagnostics_interp_vert.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
diagnostics_interp_vert.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_interp_vert.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
diagnostics_interp_vert.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
diagnostics_interp_vert.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
diagnostics_interp_vert.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagnostics_interp_vert.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
diagnostics_interp_vert.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
diagnostics_interp_vert.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diagnostics_interp_vert.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
diagnostics_interp_vert.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diagnostics_is_on.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_is_on.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagnostics_is_on.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagnostics_list_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_list_check.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_list_check.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_list_check.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_main_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_main_init.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAG_OPTIONS.h
diagnostics_main_init.f: EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_mnc_out.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_mnc_out.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_mnc_out.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h DXC_MACROS.h
diagnostics_mnc_out.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diagnostics_mnc_out.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diagnostics_mnc_out.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
diagnostics_mnc_out.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_mnc_out.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagnostics_mnc_out.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_mnc_out.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
diagnostics_mnc_out.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diagnostics_mnc_out.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_mnc_out.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_mnc_out.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_mnc_out.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagnostics_mnc_out.f: YC_MACROS.h
diagnostics_out.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_out.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagnostics_out.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diagnostics_out.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diagnostics_out.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
diagnostics_out.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_out.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagnostics_out.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagnostics_out.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagnostics_out.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagnostics_out.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagnostics_out.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagnostics_out.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagnostics_out.f: YC_MACROS.h
diagnostics_read_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_read_pickup.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_read_pickup.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_read_pickup.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_readparms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_readparms.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_readparms.f: DIAGNOSTICS_CALC.h DIAGNOSTICS_P2SHARE.h
diagnostics_readparms.f: DIAGNOSTICS_SIZE.h DIAGSTATS_REGIONS.h DIAG_OPTIONS.h
diagnostics_readparms.f: EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_scale_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_scale_fill.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_scale_fill.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_scale_fill.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_scale_fill_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_scale_fill_rs.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_scale_fill_rs.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_scale_fill_rs.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_set_calc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_set_calc.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS_CALC.h
diagnostics_set_calc.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diagnostics_set_calc.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diagnostics_set_calc.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
diagnostics_set_calc.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diagnostics_set_calc.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diagnostics_set_calc.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
diagnostics_set_calc.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
diagnostics_set_calc.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagnostics_set_calc.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
diagnostics_set_calc.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
diagnostics_set_calc.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diagnostics_set_calc.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
diagnostics_set_calc.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diagnostics_set_levels.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_set_levels.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_set_levels.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_set_levels.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_set_pointers.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_set_pointers.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_set_pointers.f: DIAGNOSTICS_P2SHARE.h DIAGNOSTICS_SIZE.h
diagnostics_set_pointers.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagnostics_setdiag.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_setdiag.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_setdiag.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_setdiag.f: PACKAGES_CONFIG.h SIZE.h
diagnostics_status_error.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_status_error.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_status_error.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_status_error.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_sum_levels.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_sum_levels.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_sum_levels.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h DXC_MACROS.h
diagnostics_sum_levels.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diagnostics_sum_levels.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diagnostics_sum_levels.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
diagnostics_sum_levels.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagnostics_sum_levels.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
diagnostics_sum_levels.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
diagnostics_sum_levels.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
diagnostics_sum_levels.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagnostics_sum_levels.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
diagnostics_sum_levels.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
diagnostics_sum_levels.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diagnostics_sum_levels.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
diagnostics_sum_levels.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diagnostics_summary.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_summary.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_summary.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_summary.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_switch_onoff.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_switch_onoff.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_switch_onoff.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_switch_onoff.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_utils.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_utils.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagnostics_utils.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
diagnostics_utils.f: SIZE.h
diagnostics_write.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagnostics_write.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_P2SHARE.h
diagnostics_write.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_write.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_write_adj.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_write_adj.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_write_adj.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_write_adj.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagnostics_write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagnostics_write_pickup.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagnostics_write_pickup.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagnostics_write_pickup.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diags_oceanic_surf_flux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diags_oceanic_surf_flux.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
diags_oceanic_surf_flux.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
diags_oceanic_surf_flux.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
diags_oceanic_surf_flux.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
diags_oceanic_surf_flux.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diags_oceanic_surf_flux.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diags_oceanic_surf_flux.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diags_oceanic_surf_flux.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diags_oceanic_surf_flux.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
diags_oceanic_surf_flux.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diags_oceanic_surf_flux.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diags_oceanic_surf_flux.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diags_oceanic_surf_flux.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diags_oceanic_surf_flux.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
diags_oceanic_surf_flux.f: XC_MACROS.h YC_MACROS.h
diags_phi_hyd.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diags_phi_hyd.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diags_phi_hyd.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
diags_phi_hyd.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
diags_phi_hyd.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diags_phi_hyd.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diags_phi_hyd.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diags_phi_hyd.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diags_phi_hyd.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diags_phi_hyd.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diags_phi_hyd.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
diags_phi_hyd.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diags_phi_rlow.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diags_phi_rlow.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diags_phi_rlow.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diags_phi_rlow.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
diags_phi_rlow.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
diags_phi_rlow.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
diags_phi_rlow.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diags_phi_rlow.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diags_phi_rlow.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diags_phi_rlow.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diags_phi_rlow.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
diags_phi_rlow.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diags_rho.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diags_rho.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
diags_rho.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
diags_rho.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diags_rho.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
diags_rho.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diags_rho.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diags_rho.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diags_rho.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diags_rho.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
diags_rho.f: XC_MACROS.h YC_MACROS.h
diags_sound_speed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diags_sound_speed.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diags_sound_speed.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diags_sound_speed.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h EOS.h FCORI_MACROS.h
diags_sound_speed.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diags_sound_speed.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diags_sound_speed.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
diags_sound_speed.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
diags_sound_speed.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diags_sound_speed.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
diags_sound_speed.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diags_sound_speed.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
diags_sound_speed.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
diagstats_ascii_out.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagstats_ascii_out.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagstats_ascii_out.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagstats_ascii_out.f: PACKAGES_CONFIG.h SIZE.h
diagstats_calc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_calc.f: CPP_OPTIONS.h DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h
diagstats_calc.f: SIZE.h
diagstats_clear.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_clear.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_clear.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagstats_close_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagstats_close_io.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagstats_close_io.f: DIAGNOSTICS_SIZE.h DIAG_OPTIONS.h EEPARAMS.h
diagstats_close_io.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagstats_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_fill.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_fill.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagstats_global.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_global.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_global.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagstats_ini_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_ini_io.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_ini_io.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
diagstats_ini_io.f: SIZE.h
diagstats_local.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_local.f: CPP_OPTIONS.h DIAGNOSTICS_SIZE.h DIAGSTATS_REGIONS.h
diagstats_local.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diagstats_local.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diagstats_local.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
diagstats_local.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
diagstats_local.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
diagstats_local.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagstats_local.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagstats_local.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagstats_local.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagstats_local.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagstats_local.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagstats_local.f: YC_MACROS.h
diagstats_mnc_out.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_mnc_out.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_mnc_out.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diagstats_mnc_out.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diagstats_mnc_out.f: DYU_MACROS.h EEPARAMS.h EESUPPORT.h FCORI_MACROS.h GRID.h
diagstats_mnc_out.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
diagstats_mnc_out.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
diagstats_mnc_out.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
diagstats_mnc_out.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
diagstats_mnc_out.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
diagstats_mnc_out.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
diagstats_mnc_out.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
diagstats_mnc_out.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
diagstats_mnc_out.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h mpif.h
diagstats_others_calc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagstats_others_calc.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAG_OPTIONS.h
diagstats_others_calc.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagstats_output.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_output.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_output.f: DIAG_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
diagstats_output.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
diagstats_output.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
diagstats_output.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
diagstats_output.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
diagstats_output.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
diagstats_output.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
diagstats_output.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
diagstats_output.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
diagstats_output.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
diagstats_output.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
diagstats_output.f: YC_MACROS.h
diagstats_set_pointers.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagstats_set_pointers.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS.h
diagstats_set_pointers.f: DIAGNOSTICS_SIZE.h DIAGSTATS_REGIONS.h
diagstats_set_pointers.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diagstats_set_regions.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diagstats_set_regions.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DIAGNOSTICS_SIZE.h
diagstats_set_regions.f: DIAGSTATS_REGIONS.h DIAG_OPTIONS.h EEPARAMS.h
diagstats_set_regions.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
diagstats_setdiag.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
diagstats_setdiag.f: CPP_OPTIONS.h DIAGNOSTICS.h DIAGNOSTICS_SIZE.h
diagstats_setdiag.f: DIAG_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
diff_phase_multiple.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
diff_phase_multiple.f: CPP_EEOPTIONS.h
different_multiple.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
different_multiple.f: CPP_EEOPTIONS.h
do_atmospheric_phys.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
do_atmospheric_phys.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
do_atmospheric_phys.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
do_atmospheric_phys.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
do_atmospheric_phys.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
do_atmospheric_phys.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
do_atmospheric_phys.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
do_atmospheric_phys.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
do_atmospheric_phys.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
do_atmospheric_phys.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
do_atmospheric_phys.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
do_atmospheric_phys.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
do_atmospheric_phys.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
do_atmospheric_phys.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
do_atmospheric_phys.f: YC_MACROS.h
do_fields_blocking_exchanges.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
do_fields_blocking_exchanges.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DYNVARS.h
do_fields_blocking_exchanges.f: EEPARAMS.h GAD.h PACKAGES_CONFIG.h PARAMS.h
do_fields_blocking_exchanges.f: SIZE.h
do_oceanic_phys.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
do_oceanic_phys.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
do_oceanic_phys.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
do_oceanic_phys.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
do_oceanic_phys.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
do_oceanic_phys.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
do_oceanic_phys.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
do_oceanic_phys.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
do_oceanic_phys.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
do_oceanic_phys.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
do_oceanic_phys.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
do_oceanic_phys.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
do_oceanic_phys.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
do_stagger_fields_exchanges.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
do_stagger_fields_exchanges.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DYNVARS.h
do_stagger_fields_exchanges.f: EEPARAMS.h PARAMS.h SIZE.h
do_statevars_diags.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
do_statevars_diags.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
do_statevars_diags.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
do_the_model_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
do_the_model_io.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h PACKAGES_CONFIG.h
do_the_model_io.f: PARAMS.h SIZE.h
do_write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
do_write_pickup.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
do_write_pickup.f: RESTART.h SIZE.h
dynamics.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
dynamics.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
dynamics.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
dynamics.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
dynamics.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
dynamics.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_VISC.h PACKAGES_CONFIG.h
dynamics.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
dynamics.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
dynamics.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
dynamics.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
dynamics.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
dynamics.f: XC_MACROS.h YC_MACROS.h
eeboot.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h EEPARAMS.h
eeboot.f: EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
eeboot_minimal.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
eeboot_minimal.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
eedata_example.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
eedata_example.f: EEPARAMS.h SIZE.h
eedie.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h EEPARAMS.h
eedie.f: EESUPPORT.h SIZE.h mpif.h
eeintro_msg.f: /usr/include/stdc-predef.h BUILD_INFO.h CPP_EEMACROS.h
eeintro_msg.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
eeset_parms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
eeset_parms.f: EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h mpif.h
eewrite_eeenv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
eewrite_eeenv.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
exch0_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch0_r4.f: EEPARAMS.h SIZE.h
exch0_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch0_r8.f: EEPARAMS.h SIZE.h
exch0_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch0_rl.f: EEPARAMS.h SIZE.h
exch0_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch0_rs.f: EEPARAMS.h SIZE.h
exch1_bg_r4_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_bg_r4_cube.f: EEPARAMS.h SIZE.h
exch1_bg_r8_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_bg_r8_cube.f: EEPARAMS.h SIZE.h
exch1_bg_rl_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_bg_rl_cube.f: EEPARAMS.h SIZE.h
exch1_bg_rs_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_bg_rs_cube.f: EEPARAMS.h SIZE.h
exch1_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_r4.f: EEPARAMS.h EXCH.h SIZE.h
exch1_r4_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_r4_cube.f: EEPARAMS.h SIZE.h
exch1_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_r8.f: EEPARAMS.h EXCH.h SIZE.h
exch1_r8_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_r8_cube.f: EEPARAMS.h SIZE.h
exch1_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rl.f: EEPARAMS.h EXCH.h SIZE.h
exch1_rl_ad.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rl_ad.f: EEPARAMS.h EXCH.h SIZE.h
exch1_rl_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rl_cube.f: EEPARAMS.h SIZE.h
exch1_rl_cube_ad.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rl_cube_ad.f: EEPARAMS.h SIZE.h
exch1_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rs.f: EEPARAMS.h EXCH.h SIZE.h
exch1_rs_ad.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rs_ad.f: EEPARAMS.h EXCH.h SIZE.h
exch1_rs_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rs_cube.f: EEPARAMS.h SIZE.h
exch1_rs_cube_ad.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_rs_cube_ad.f: EEPARAMS.h SIZE.h
exch1_uv_r4_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_uv_r4_cube.f: EEPARAMS.h SIZE.h
exch1_uv_r8_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_uv_r8_cube.f: EEPARAMS.h SIZE.h
exch1_uv_rl_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_uv_rl_cube.f: EEPARAMS.h SIZE.h
exch1_uv_rs_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_uv_rs_cube.f: EEPARAMS.h SIZE.h
exch1_z_r4_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_z_r4_cube.f: EEPARAMS.h SIZE.h
exch1_z_r8_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_z_r8_cube.f: EEPARAMS.h SIZE.h
exch1_z_rl_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_z_rl_cube.f: EEPARAMS.h SIZE.h
exch1_z_rs_cube.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch1_z_rs_cube.f: EEPARAMS.h SIZE.h
exch_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_3d_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_3d_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_3d_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_3d_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_cycle_ebl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_cycle_ebl.f: EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h mpif.h
exch_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_init.f: EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h mpif.h
exch_r4_recv_get_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r4_recv_get_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r4_recv_get_x.f: mpif.h
exch_r4_recv_get_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r4_recv_get_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r4_recv_get_y.f: mpif.h
exch_r4_send_put_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r4_send_put_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r4_send_put_x.f: mpif.h
exch_r4_send_put_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r4_send_put_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r4_send_put_y.f: mpif.h
exch_r8_recv_get_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r8_recv_get_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r8_recv_get_x.f: mpif.h
exch_r8_recv_get_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r8_recv_get_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r8_recv_get_y.f: mpif.h
exch_r8_send_put_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r8_send_put_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r8_send_put_x.f: mpif.h
exch_r8_send_put_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_r8_send_put_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_r8_send_put_y.f: mpif.h
exch_rl_recv_get_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rl_recv_get_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rl_recv_get_x.f: mpif.h
exch_rl_recv_get_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rl_recv_get_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rl_recv_get_y.f: mpif.h
exch_rl_send_put_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rl_send_put_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rl_send_put_x.f: mpif.h
exch_rl_send_put_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rl_send_put_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rl_send_put_y.f: mpif.h
exch_rs_recv_get_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rs_recv_get_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rs_recv_get_x.f: mpif.h
exch_rs_recv_get_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rs_recv_get_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rs_recv_get_y.f: mpif.h
exch_rs_send_put_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rs_send_put_x.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rs_send_put_x.f: mpif.h
exch_rs_send_put_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_rs_send_put_y.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h SIZE.h
exch_rs_send_put_y.f: mpif.h
exch_s3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_s3d_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_s3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_s3d_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_s3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_s3d_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_s3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_s3d_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_sm_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_sm_3d_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_sm_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_sm_3d_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_sm_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_sm_3d_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_sm_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_sm_3d_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_tap_b.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_tap_b.f: EEPARAMS.h EXCH.h PACKAGES_CONFIG.h SIZE.h
exch_tap_d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_tap_d.f: EEPARAMS.h EXCH.h PACKAGES_CONFIG.h SIZE.h
exch_uv_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_3d_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_3d_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_3d_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_3d_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_agrid_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_agrid_3d_r4.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_agrid_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_agrid_3d_r8.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_agrid_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_agrid_3d_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_agrid_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_agrid_3d_rs.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_bgrid_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_bgrid_3d_r4.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_bgrid_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_bgrid_3d_r8.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_bgrid_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_bgrid_3d_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_bgrid_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_bgrid_3d_rs.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_dgrid_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_dgrid_3d_r4.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_dgrid_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_dgrid_3d_r8.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_dgrid_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_dgrid_3d_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_dgrid_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
exch_uv_dgrid_3d_rs.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xy_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xy_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xy_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xy_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xy_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xy_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xy_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xy_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xyz_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xyz_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xyz_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xyz_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xyz_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xyz_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_uv_xyz_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_uv_xyz_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xy_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xy_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xy_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xy_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xy_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xy_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xy_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xy_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xyz_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xyz_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xyz_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xyz_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xyz_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xyz_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_xyz_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_xyz_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_z_3d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_z_3d_r4.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_z_3d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_z_3d_r8.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_z_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_z_3d_rl.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
exch_z_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
exch_z_3d_rs.f: EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
external_fields_load.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
external_fields_load.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
external_fields_load.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
external_fields_load.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
external_fields_load.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
external_fields_load.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
external_fields_load.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
external_fields_load.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
external_fields_load.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
external_fields_load.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
external_fields_load.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
external_fields_load.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
external_fields_load.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
external_fields_load.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
external_fields_load.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
external_forcing.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
external_forcing.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
external_forcing.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
external_forcing.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
external_forcing.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
external_forcing.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
external_forcing.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
external_forcing.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
external_forcing.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
external_forcing.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
external_forcing.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
external_forcing.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
external_forcing.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
external_forcing_surf.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
external_forcing_surf.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
external_forcing_surf.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
external_forcing_surf.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
external_forcing_surf.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
external_forcing_surf.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
external_forcing_surf.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
external_forcing_surf.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
external_forcing_surf.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
external_forcing_surf.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
external_forcing_surf.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
external_forcing_surf.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
external_forcing_surf.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
external_forcing_surf.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
external_forcing_surf.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
external_forcing_surf.f: YC_MACROS.h
fill_cs_corner_ag_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
fill_cs_corner_ag_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
fill_cs_corner_tr_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
fill_cs_corner_tr_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
fill_cs_corner_uv_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
fill_cs_corner_uv_rl.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
fill_cs_corner_uv_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
fill_cs_corner_uv_rs.f: CPP_EEOPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
fill_in_corners_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
fill_in_corners_rl.f: CPP_EEOPTIONS.h EEPARAMS.h SIZE.h
find_alpha.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
find_alpha.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
find_alpha.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
find_alpha.f: DYU_MACROS.h EEPARAMS.h EOS.h FCORI_MACROS.h GRID.h
find_alpha.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
find_alpha.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
find_alpha.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
find_alpha.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
find_alpha.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
find_alpha.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
find_alpha.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
find_alpha.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
find_hyd_press_1d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
find_hyd_press_1d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
find_hyd_press_1d.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
find_hyd_press_1d.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
find_hyd_press_1d.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
find_hyd_press_1d.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
find_hyd_press_1d.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
find_hyd_press_1d.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
find_hyd_press_1d.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
find_hyd_press_1d.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
find_hyd_press_1d.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
find_hyd_press_1d.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
find_hyd_press_1d.f: YC_MACROS.h
find_rho.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
find_rho.f: CPP_OPTIONS.h EEPARAMS.h EOS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
fool_the_compiler.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
forcing_surf_relax.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
forcing_surf_relax.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
forcing_surf_relax.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
forcing_surf_relax.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
forcing_surf_relax.f: FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
forcing_surf_relax.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
forcing_surf_relax.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
forcing_surf_relax.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
forcing_surf_relax.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
forcing_surf_relax.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
forcing_surf_relax.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
forcing_surf_relax.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
forcing_surf_relax.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
forcing_surf_relax.f: XC_MACROS.h YC_MACROS.h
forward_step.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
forward_step.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h GAD_OPTIONS.h
forward_step.f: PACKAGES_CONFIG.h PARAMS.h SIGREG.h SIZE.h
freesurf_rescale_g.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
freesurf_rescale_g.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
freesurf_rescale_g.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
freesurf_rescale_g.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
freesurf_rescale_g.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
freesurf_rescale_g.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
freesurf_rescale_g.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
freesurf_rescale_g.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
freesurf_rescale_g.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
freesurf_rescale_g.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
freesurf_rescale_g.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
freesurf_rescale_g.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
freesurf_rescale_g.f: XC_MACROS.h YC_MACROS.h
freeze_surface.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
freeze_surface.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
freeze_surface.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
freeze_surface.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
freeze_surface.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
freeze_surface.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
freeze_surface.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
freeze_surface.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
freeze_surface.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
freeze_surface.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
freeze_surface.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
freeze_surface.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
freeze_surface.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_advection.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_advection.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_advection.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_advection.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_advection.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_advection.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_advection.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_advection.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_advection.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_advection.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_advection.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_advection.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_advection.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_advscheme.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_advscheme.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_advscheme.f: PACKAGES_CONFIG.h SIZE.h
gad_biharm_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_biharm_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_biharm_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_biharm_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h
gad_biharm_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_biharm_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
gad_biharm_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_biharm_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_biharm_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_biharm_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_biharm_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_biharm_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_biharm_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_biharm_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_biharm_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_biharm_x.f: DYU_MACROS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_biharm_x.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_biharm_x.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_biharm_x.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_biharm_x.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_biharm_x.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_biharm_x.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_biharm_x.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_biharm_x.f: YC_MACROS.h
gad_biharm_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_biharm_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_biharm_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_biharm_y.f: DYU_MACROS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_biharm_y.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_biharm_y.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_biharm_y.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_biharm_y.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_biharm_y.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_biharm_y.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_biharm_y.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_biharm_y.f: YC_MACROS.h
gad_c2_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c2_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c2_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c2_adv_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h
gad_c2_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_c2_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
gad_c2_adv_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_c2_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_c2_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_c2_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_c2_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_c2_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_c2_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c2_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c2_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c2_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_c2_adv_x.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_c2_adv_x.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_c2_adv_x.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_c2_adv_x.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_c2_adv_x.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_c2_adv_x.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_c2_adv_x.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_c2_adv_x.f: YC_MACROS.h
gad_c2_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c2_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c2_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c2_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_c2_adv_y.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_c2_adv_y.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_c2_adv_y.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_c2_adv_y.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_c2_adv_y.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_c2_adv_y.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_c2_adv_y.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_c2_adv_y.f: YC_MACROS.h
gad_c2_impl_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c2_impl_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c2_impl_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c2_impl_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h
gad_c2_impl_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_c2_impl_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
gad_c2_impl_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_c2_impl_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_c2_impl_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_c2_impl_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_c2_impl_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_c2_impl_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_c4_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c4_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c4_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c4_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_c4_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_c4_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_c4_adv_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_c4_adv_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_c4_adv_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_c4_adv_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_c4_adv_r.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_c4_adv_r.f: YC_MACROS.h
gad_c4_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c4_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c4_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c4_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_c4_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_c4_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_c4_adv_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_c4_adv_x.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_c4_adv_x.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_c4_adv_x.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_c4_adv_x.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_c4_adv_x.f: YC_MACROS.h
gad_c4_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_c4_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_c4_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_c4_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_c4_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_c4_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_c4_adv_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_c4_adv_y.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_c4_adv_y.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_c4_adv_y.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_c4_adv_y.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_c4_adv_y.f: YC_MACROS.h
gad_calc_rhs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_calc_rhs.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_calc_rhs.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_calc_rhs.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_calc_rhs.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_calc_rhs.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_calc_rhs.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_calc_rhs.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
gad_calc_rhs.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_calc_rhs.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_calc_rhs.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
gad_calc_rhs.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_check.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h PACKAGES_CONFIG.h
gad_check.f: PARAMS.h SIZE.h
gad_del2.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_del2.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
gad_del2.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h FCORI_MACROS.h
gad_del2.f: GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_del2.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_del2.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_del2.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_del2.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_del2.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_del2.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
gad_del2.f: XC_MACROS.h YC_MACROS.h
gad_diagnostics_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_diagnostics_init.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h GAD.h
gad_diagnostics_init.f: GAD_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
gad_diagnostics_state.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_diagnostics_state.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
gad_diagnostics_state.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
gad_diagnostics_state.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
gad_diagnostics_state.f: FCORI_MACROS.h GAD.h GAD_OPTIONS.h GAD_SOM_VARS.h
gad_diagnostics_state.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_diagnostics_state.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_diagnostics_state.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_diagnostics_state.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_diagnostics_state.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_diagnostics_state.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_diagnostics_state.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_diagnostics_state.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_diagnostics_state.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_diagnostics_state.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_diff_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_diff_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_diff_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
gad_diff_r.f: EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_diff_r.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_diff_r.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
gad_diff_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_diff_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_diff_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_diff_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_diff_r.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_diff_r.f: YC_MACROS.h
gad_diff_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_diff_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_diff_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
gad_diff_x.f: FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
gad_diff_x.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_diff_x.f: PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_diff_x.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
gad_diff_x.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_diff_x.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_diff_x.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_diff_x.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_diff_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_diff_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_diff_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
gad_diff_y.f: FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
gad_diff_y.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_diff_y.f: PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_diff_y.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
gad_diff_y.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_diff_y.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_diff_y.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_diff_y.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_dst2u1_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst2u1_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst2u1_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst2u1_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst2u1_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst2u1_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst2u1_adv_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst2u1_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst2u1_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst2u1_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst2u1_adv_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst2u1_adv_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst2u1_adv_r.f: YC_MACROS.h
gad_dst2u1_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst2u1_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst2u1_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst2u1_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst2u1_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst2u1_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst2u1_adv_x.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst2u1_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst2u1_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst2u1_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst2u1_adv_x.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst2u1_adv_x.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst2u1_adv_x.f: YC_MACROS.h
gad_dst2u1_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst2u1_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst2u1_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst2u1_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst2u1_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst2u1_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst2u1_adv_y.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst2u1_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst2u1_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst2u1_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst2u1_adv_y.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst2u1_adv_y.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst2u1_adv_y.f: YC_MACROS.h
gad_dst2u1_impl_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst2u1_impl_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst2u1_impl_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst2u1_impl_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h
gad_dst2u1_impl_r.f: GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
gad_dst2u1_impl_r.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_dst2u1_impl_r.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
gad_dst2u1_impl_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst2u1_impl_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst2u1_impl_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst2u1_impl_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst2u1_impl_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst2u1_impl_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst2u1_impl_r.f: YC_MACROS.h
gad_dst3_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_dst3_adv_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_dst3_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_dst3_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_dst3_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_dst3_adv_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_dst3_adv_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_dst3_adv_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_dst3_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_dst3_adv_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_dst3_adv_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_dst3_adv_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_dst3fl_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3fl_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3fl_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3fl_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3fl_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3fl_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst3fl_adv_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3fl_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3fl_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3fl_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst3fl_adv_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst3fl_adv_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst3fl_adv_r.f: YC_MACROS.h
gad_dst3fl_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3fl_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3fl_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3fl_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3fl_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3fl_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst3fl_adv_x.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3fl_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3fl_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3fl_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst3fl_adv_x.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst3fl_adv_x.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst3fl_adv_x.f: YC_MACROS.h
gad_dst3fl_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3fl_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3fl_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3fl_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_dst3fl_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_dst3fl_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_dst3fl_adv_y.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3fl_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3fl_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3fl_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst3fl_adv_y.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst3fl_adv_y.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst3fl_adv_y.f: YC_MACROS.h
gad_dst3fl_impl_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_dst3fl_impl_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_dst3fl_impl_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_dst3fl_impl_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h
gad_dst3fl_impl_r.f: GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
gad_dst3fl_impl_r.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_dst3fl_impl_r.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
gad_dst3fl_impl_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_dst3fl_impl_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_dst3fl_impl_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_dst3fl_impl_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_dst3fl_impl_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_dst3fl_impl_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_dst3fl_impl_r.f: YC_MACROS.h
gad_exch_som.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_exch_som.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h PACKAGES_CONFIG.h
gad_exch_som.f: SIZE.h
gad_fluxlimit_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_fluxlimit_adv_r.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
gad_fluxlimit_adv_r.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
gad_fluxlimit_adv_r.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
gad_fluxlimit_adv_r.f: GAD_FLUX_LIMITER.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_fluxlimit_adv_r.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_fluxlimit_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_fluxlimit_adv_r.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_fluxlimit_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_fluxlimit_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_fluxlimit_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_fluxlimit_adv_x.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
gad_fluxlimit_adv_x.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
gad_fluxlimit_adv_x.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
gad_fluxlimit_adv_x.f: GAD_FLUX_LIMITER.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_fluxlimit_adv_x.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_fluxlimit_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_fluxlimit_adv_x.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_fluxlimit_adv_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_fluxlimit_adv_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_fluxlimit_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_fluxlimit_adv_y.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
gad_fluxlimit_adv_y.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
gad_fluxlimit_adv_y.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
gad_fluxlimit_adv_y.f: GAD_FLUX_LIMITER.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_fluxlimit_adv_y.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_fluxlimit_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_fluxlimit_adv_y.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_fluxlimit_adv_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_fluxlimit_adv_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_fluxlimit_impl_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_fluxlimit_impl_r.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
gad_fluxlimit_impl_r.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
gad_fluxlimit_impl_r.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
gad_fluxlimit_impl_r.f: FCORI_MACROS.h GAD_FLUX_LIMITER.h GAD_OPTIONS.h GRID.h
gad_fluxlimit_impl_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_fluxlimit_impl_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_fluxlimit_impl_r.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_fluxlimit_impl_r.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_fluxlimit_impl_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_fluxlimit_impl_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
gad_fluxlimit_impl_r.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_fluxlimit_impl_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_fluxlimit_impl_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_fluxlimit_impl_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_grad_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_grad_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_grad_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
gad_grad_x.f: EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_grad_x.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_grad_x.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_grad_x.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_grad_x.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_grad_x.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_grad_x.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_grad_x.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_grad_x.f: YC_MACROS.h
gad_grad_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_grad_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_grad_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
gad_grad_y.f: EEPARAMS.h FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h
gad_grad_y.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
gad_grad_y.f: MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h
gad_grad_y.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_grad_y.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_grad_y.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_grad_y.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_grad_y.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_grad_y.f: YC_MACROS.h
gad_implicit_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_implicit_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_implicit_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_implicit_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_implicit_r.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_implicit_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_implicit_r.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_implicit_r.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_implicit_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_implicit_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_implicit_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_implicit_r.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
gad_implicit_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_init_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_init_fixed.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_init_fixed.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
gad_init_varia.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_init_varia.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h GAD_SOM_VARS.h
gad_init_varia.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
gad_os7mp_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_os7mp_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_os7mp_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_os7mp_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_os7mp_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_os7mp_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_os7mp_adv_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_os7mp_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_os7mp_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_os7mp_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_os7mp_adv_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_os7mp_adv_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_os7mp_adv_r.f: YC_MACROS.h
gad_os7mp_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_os7mp_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_os7mp_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_os7mp_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_os7mp_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_os7mp_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_os7mp_adv_x.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_os7mp_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_os7mp_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_os7mp_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_os7mp_adv_x.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_os7mp_adv_x.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_os7mp_adv_x.f: YC_MACROS.h
gad_os7mp_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_os7mp_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_os7mp_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_os7mp_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_os7mp_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_os7mp_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
gad_os7mp_adv_y.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_os7mp_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_os7mp_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_os7mp_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
gad_os7mp_adv_y.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
gad_os7mp_adv_y.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_os7mp_adv_y.f: YC_MACROS.h
gad_osc_hat_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_hat_r.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_osc_hat_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_hat_x.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_osc_hat_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_hat_y.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_osc_mul_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_mul_r.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_osc_mul_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_mul_x.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_osc_mul_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_osc_mul_y.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_plm_fun.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_plm_fun.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h
gad_ppm_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_adv_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_adv_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_adv_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_adv_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_adv_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_adv_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_adv_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_flx_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_flx_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_flx_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_flx_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_flx_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_flx_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_flx_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_flx_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_flx_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_flx_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_flx_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_flx_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_flx_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_flx_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_flx_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_flx_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_flx_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_flx_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_flx_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_flx_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_flx_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_flx_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_flx_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_flx_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_flx_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_flx_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_flx_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_flx_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_flx_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_flx_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_flx_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_flx_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_flx_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_flx_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_flx_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_flx_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_fun.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_fun.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h
gad_ppm_hat_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_hat_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_hat_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_hat_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_hat_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_hat_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_hat_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_hat_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_hat_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_hat_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_hat_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_hat_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_hat_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_hat_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_hat_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_hat_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_hat_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_hat_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_hat_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_hat_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_hat_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_hat_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_hat_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_hat_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_hat_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_hat_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_hat_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_hat_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_hat_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_hat_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_hat_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_hat_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_hat_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_hat_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_hat_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_hat_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_p3e_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_p3e_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_p3e_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_p3e_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_p3e_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_p3e_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_p3e_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_p3e_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_p3e_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_p3e_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_p3e_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_p3e_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_p3e_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_p3e_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_p3e_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_p3e_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_p3e_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_p3e_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_p3e_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_p3e_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_p3e_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_p3e_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_p3e_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_p3e_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_ppm_p3e_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_ppm_p3e_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_ppm_p3e_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_ppm_p3e_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_ppm_p3e_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_ppm_p3e_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_ppm_p3e_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_ppm_p3e_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_ppm_p3e_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_ppm_p3e_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_ppm_p3e_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_ppm_p3e_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_adv_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_adv_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_adv_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_adv_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_adv_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_adv_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_adv_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_adv_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_adv_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_adv_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_adv_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_adv_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_adv_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_adv_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_adv_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_adv_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_adv_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_adv_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_adv_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_adv_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_adv_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_adv_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_adv_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_adv_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_adv_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_adv_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_flx_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_flx_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_flx_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_flx_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_flx_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_flx_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_flx_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_flx_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_flx_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_flx_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_flx_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_flx_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_flx_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_flx_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_flx_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_flx_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_flx_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_flx_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_flx_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_flx_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_flx_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_flx_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_flx_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_flx_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_flx_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_flx_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_flx_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_flx_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_flx_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_flx_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_flx_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_flx_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_flx_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_flx_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_flx_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_flx_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_fun.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_fun.f: CPP_OPTIONS.h GAD_OPTIONS.h PACKAGES_CONFIG.h
gad_pqm_hat_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_hat_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_hat_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_hat_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_hat_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_hat_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_hat_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_hat_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_hat_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_hat_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_hat_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_hat_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_hat_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_hat_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_hat_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_hat_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_hat_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_hat_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_hat_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_hat_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_hat_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_hat_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_hat_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_hat_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_hat_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_hat_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_hat_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_hat_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_hat_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_hat_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_hat_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_hat_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_hat_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_hat_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_hat_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_hat_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_p5e_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_p5e_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_p5e_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_p5e_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_p5e_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_p5e_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_p5e_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_p5e_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_p5e_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_p5e_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_p5e_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_p5e_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_p5e_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_p5e_x.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_p5e_x.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_p5e_x.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_p5e_x.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_p5e_x.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_p5e_x.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_p5e_x.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_p5e_x.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_p5e_x.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_p5e_x.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_p5e_x.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_pqm_p5e_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_pqm_p5e_y.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_pqm_p5e_y.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_pqm_p5e_y.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_pqm_p5e_y.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_pqm_p5e_y.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_pqm_p5e_y.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
gad_pqm_p5e_y.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
gad_pqm_p5e_y.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
gad_pqm_p5e_y.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_pqm_p5e_y.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_pqm_p5e_y.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_read_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_read_pickup.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h GAD_SOM_VARS.h
gad_read_pickup.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
gad_som_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_som_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_som_adv_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_som_adv_r.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_som_adv_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_som_adv_r.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_som_adv_r.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_som_adv_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_som_adv_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_som_adv_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_som_adv_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_som_adv_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_som_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_adv_x.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_som_adv_x.f: PACKAGES_CONFIG.h SIZE.h
gad_som_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_adv_y.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_som_adv_y.f: PACKAGES_CONFIG.h SIZE.h
gad_som_advect.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_advect.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_som_advect.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_som_advect.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_som_advect.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_som_advect.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_som_advect.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_som_advect.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_som_advect.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_som_advect.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_som_advect.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_som_advect.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_som_advect.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_som_exchanges.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_exchanges.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_som_exchanges.f: GAD_SOM_VARS.h PACKAGES_CONFIG.h SIZE.h
gad_som_fill_cs_corner.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_som_fill_cs_corner.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h GAD.h
gad_som_fill_cs_corner.f: GAD_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
gad_som_lim_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_som_lim_r.f: CPP_OPTIONS.h GAD.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_som_prep_cs_corner.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
gad_som_prep_cs_corner.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h GAD.h
gad_som_prep_cs_corner.f: GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_u3_adv_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_u3_adv_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_u3_adv_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_u3_adv_r.f: DYU_MACROS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h GRID.h
gad_u3_adv_r.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
gad_u3_adv_r.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h RAS_MACROS.h
gad_u3_adv_r.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_u3_adv_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_u3_adv_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_u3_adv_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
gad_u3_adv_r.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
gad_u3_adv_r.f: YC_MACROS.h
gad_u3_adv_x.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_u3_adv_x.f: CPP_OPTIONS.h GAD.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_u3_adv_y.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_u3_adv_y.f: CPP_OPTIONS.h GAD.h GAD_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
gad_u3c4_impl_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_u3c4_impl_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
gad_u3c4_impl_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
gad_u3c4_impl_r.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h GAD_OPTIONS.h
gad_u3c4_impl_r.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
gad_u3c4_impl_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
gad_u3c4_impl_r.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
gad_u3c4_impl_r.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
gad_u3c4_impl_r.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
gad_u3c4_impl_r.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
gad_u3c4_impl_r.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
gad_u3c4_impl_r.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
gad_u3c4_impl_r.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gad_write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gad_write_pickup.f: CPP_OPTIONS.h EEPARAMS.h GAD.h GAD_OPTIONS.h
gad_write_pickup.f: GAD_SOM_VARS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
gather_2d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_2d_r4.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
gather_2d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_2d_r8.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
gather_2d_wh_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_2d_wh_r4.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
gather_2d_wh_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_2d_wh_r8.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
gather_vec_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_vec_r4.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
gather_vec_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_vec_r8.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
gather_xz.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_xz.f: CPP_OPTIONS.h EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
gather_yz.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gather_yz.f: CPP_OPTIONS.h EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
get_periodic_interval.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
get_periodic_interval.f: CPP_EEOPTIONS.h EEPARAMS.h
get_write_global_fld.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
get_write_global_fld.f: CPP_EEOPTIONS.h
global_max.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
global_max.f: EEPARAMS.h EESUPPORT.h GLOBAL_MAX.h SIZE.h mpif.h
global_sum.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
global_sum.f: EEPARAMS.h EESUPPORT.h GLOBAL_SUM.h SIZE.h mpif.h
global_sum_singlecpu.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
global_sum_singlecpu.f: CPP_EEOPTIONS.h EEBUFF_SCPU.h EEPARAMS.h EESUPPORT.h
global_sum_singlecpu.f: GLOBAL_SUM.h PACKAGES_CONFIG.h SIZE.h mpif.h
global_sum_tile.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
global_sum_tile.f: EEPARAMS.h EESUPPORT.h GLOBAL_SUM.h SIZE.h mpif.h
global_sum_vector.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
global_sum_vector.f: EEPARAMS.h EESUPPORT.h GLOBAL_SUM.h SIZE.h mpif.h
grad_sigma.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
grad_sigma.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
grad_sigma.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
grad_sigma.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
grad_sigma.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
grad_sigma.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
grad_sigma.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
grad_sigma.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
grad_sigma.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
grad_sigma.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
grad_sigma.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
gsw_teos10.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
gsw_teos10.f: CPP_OPTIONS.h
impldiff.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
impldiff.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
impldiff.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
impldiff.f: EEPARAMS.h FCORI_MACROS.h GAD.h GRID.h GRID_MACROS.h
impldiff.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
impldiff.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
impldiff.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
impldiff.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
impldiff.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
impldiff.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
impldiff.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
impldiff.f: YC_MACROS.h
ini_cartesian_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_cartesian_grid.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
ini_cartesian_grid.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
ini_cartesian_grid.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
ini_cartesian_grid.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_cartesian_grid.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
ini_cartesian_grid.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_cartesian_grid.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_cartesian_grid.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_cartesian_grid.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_cartesian_grid.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_cartesian_grid.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
ini_cartesian_grid.f: YC_MACROS.h
ini_cg2d.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_cg2d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_cg2d.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_cg2d.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_cg2d.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
ini_cg2d.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_cg2d.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_cg2d.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_cg2d.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_cg2d.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
ini_cg2d.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_cg3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_cg3d.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
ini_communication_patterns.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_communication_patterns.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h EXCH.h
ini_communication_patterns.f: PACKAGES_CONFIG.h SIZE.h mpif.h
ini_cori.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_cori.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_cori.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_cori.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_cori.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h MONITOR.h
ini_cori.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_cori.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_cori.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_cori.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_cori.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_cori.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_curvilinear_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_curvilinear_grid.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
ini_curvilinear_grid.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
ini_curvilinear_grid.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_curvilinear_grid.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_curvilinear_grid.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_curvilinear_grid.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
ini_curvilinear_grid.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_curvilinear_grid.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_curvilinear_grid.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_curvilinear_grid.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_curvilinear_grid.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_curvilinear_grid.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_curvilinear_grid.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
ini_curvilinear_grid.f: YC_MACROS.h
ini_cylinder_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_cylinder_grid.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_cylinder_grid.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_cylinder_grid.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
ini_cylinder_grid.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_cylinder_grid.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
ini_cylinder_grid.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_cylinder_grid.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_cylinder_grid.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_cylinder_grid.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_cylinder_grid.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_cylinder_grid.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
ini_cylinder_grid.f: YC_MACROS.h
ini_depths.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_depths.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_depths.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h
ini_depths.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_depths.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_depths.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_depths.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_depths.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_depths.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_depths.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
ini_depths.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_dynvars.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_dynvars.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h MOM_COMMON_OPTIONS.h
ini_dynvars.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
ini_eos.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_eos.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_eos.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_eos.f: EOS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_eos.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_eos.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_eos.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_eos.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_eos.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_eos.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
ini_eos.f: XC_MACROS.h YC_MACROS.h
ini_ffields.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_ffields.f: CPP_OPTIONS.h EEPARAMS.h FFIELDS.h PARAMS.h SIZE.h
ini_fields.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_fields.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
ini_forcing.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_forcing.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_forcing.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_forcing.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
ini_forcing.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
ini_forcing.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
ini_forcing.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_forcing.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_forcing.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_forcing.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_forcing.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
ini_forcing.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_global_domain.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_global_domain.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_global_domain.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_global_domain.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
ini_global_domain.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_global_domain.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_global_domain.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
ini_global_domain.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_global_domain.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_global_domain.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_global_domain.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_global_domain.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
ini_global_domain.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_grid.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_grid.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_grid.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_grid.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h MONITOR.h
ini_grid.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_grid.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_grid.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_grid.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_grid.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_grid.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_linear_phisurf.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_linear_phisurf.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
ini_linear_phisurf.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
ini_linear_phisurf.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
ini_linear_phisurf.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_linear_phisurf.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_linear_phisurf.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
ini_linear_phisurf.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_linear_phisurf.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_linear_phisurf.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_linear_phisurf.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_linear_phisurf.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
ini_linear_phisurf.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_local_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_local_grid.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
ini_local_grid.f: SET_GRID.h SIZE.h
ini_masks_etc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_masks_etc.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_masks_etc.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_masks_etc.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_masks_etc.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_masks_etc.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
ini_masks_etc.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_masks_etc.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_masks_etc.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_masks_etc.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_masks_etc.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
ini_masks_etc.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_mixing.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_mixing.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_mixing.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
ini_mixing.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_mixing.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_mixing.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_mixing.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_mixing.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_mixing.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_mixing.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_mixing.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_mnc_vars.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_mnc_vars.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_mnc_vars.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_mnc_vars.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_mnc_vars.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_mnc_vars.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
ini_mnc_vars.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_mnc_vars.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_mnc_vars.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_mnc_vars.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
ini_mnc_vars.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
ini_mnc_vars.f: YC_MACROS.h
ini_model_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_model_io.f: CPP_OPTIONS.h EEBUFF_SCPU.h EEPARAMS.h EESUPPORT.h
ini_model_io.f: PACKAGES_CONFIG.h PARAMS.h RESTART.h SIZE.h mpif.h
ini_nh_fields.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_nh_fields.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_nh_fields.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_nh_fields.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_nh_fields.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_nh_fields.f: MASKW_MACROS.h NH_VARS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
ini_nh_fields.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_nh_fields.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_nh_fields.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_nh_fields.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_nh_fields.f: RECIP_HFACW_MACROS.h RESTART.h SIZE.h TANPHIATU_MACROS.h
ini_nh_fields.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_nh_vars.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_nh_vars.f: CPP_OPTIONS.h EEPARAMS.h NH_VARS.h PARAMS.h SIZE.h
ini_nlfs_vars.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_nlfs_vars.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_nlfs_vars.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_nlfs_vars.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_nlfs_vars.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_nlfs_vars.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
ini_nlfs_vars.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_nlfs_vars.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_nlfs_vars.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_nlfs_vars.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_nlfs_vars.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
ini_nlfs_vars.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_p_ground.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_p_ground.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_p_ground.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_p_ground.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_p_ground.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_p_ground.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_p_ground.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_p_ground.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_p_ground.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_p_ground.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
ini_p_ground.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_parms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_parms.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_parms.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
ini_parms.f: EOS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_parms.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_parms.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_parms.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_parms.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_parms.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_parms.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SET_GRID.h SIZE.h
ini_parms.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_pressure.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_pressure.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_pressure.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
ini_pressure.f: DYU_MACROS.h EEPARAMS.h EOS.h FCORI_MACROS.h GRID.h
ini_pressure.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
ini_pressure.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
ini_pressure.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_pressure.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_pressure.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_pressure.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_pressure.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
ini_pressure.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_procs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_procs.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
ini_psurf.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_psurf.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_psurf.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
ini_psurf.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_psurf.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_psurf.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_psurf.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_psurf.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_psurf.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_psurf.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
ini_psurf.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_salt.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_salt.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_salt.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
ini_salt.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_salt.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_salt.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_salt.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_salt.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_salt.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_salt.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_salt.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_sigma_hfac.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_sigma_hfac.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_sigma_hfac.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_sigma_hfac.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
ini_sigma_hfac.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
ini_sigma_hfac.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
ini_sigma_hfac.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_sigma_hfac.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
ini_sigma_hfac.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_sigma_hfac.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_sigma_hfac.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
ini_sigma_hfac.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_spherical_polar_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_spherical_polar_grid.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
ini_spherical_polar_grid.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_spherical_polar_grid.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_spherical_polar_grid.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
ini_spherical_polar_grid.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_spherical_polar_grid.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_spherical_polar_grid.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_spherical_polar_grid.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
ini_spherical_polar_grid.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_spherical_polar_grid.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_spherical_polar_grid.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
ini_spherical_polar_grid.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_spherical_polar_grid.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
ini_spherical_polar_grid.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_theta.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_theta.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_theta.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
ini_theta.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_theta.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_theta.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
ini_theta.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
ini_theta.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
ini_theta.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_theta.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_theta.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
ini_threading_environment.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
ini_threading_environment.f: CPP_EEOPTIONS.h EEPARAMS.h EESUPPORT.h
ini_threading_environment.f: PACKAGES_CONFIG.h SIZE.h mpif.h
ini_vel.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_vel.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
ini_vel.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
ini_vel.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
ini_vel.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
ini_vel.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_vel.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_vel.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_vel.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
ini_vel.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
ini_vel.f: XC_MACROS.h YC_MACROS.h
ini_vertical_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
ini_vertical_grid.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
ini_vertical_grid.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
ini_vertical_grid.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
ini_vertical_grid.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
ini_vertical_grid.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h PARAMS.h
ini_vertical_grid.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
ini_vertical_grid.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
ini_vertical_grid.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
ini_vertical_grid.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
ini_vertical_grid.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
ini_vertical_grid.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
ini_vertical_grid.f: YC_MACROS.h
initialise_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
initialise_fixed.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
initialise_varia.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
initialise_varia.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h PACKAGES_CONFIG.h
initialise_varia.f: PARAMS.h SIZE.h SURFACE.h
integr_continuity.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
integr_continuity.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
integr_continuity.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
integr_continuity.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
integr_continuity.f: FFIELDS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
integr_continuity.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
integr_continuity.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
integr_continuity.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
integr_continuity.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
integr_continuity.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
integr_continuity.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
integr_continuity.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
integr_continuity.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
integr_continuity.f: XC_MACROS.h YC_MACROS.h
integrate_for_w.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
integrate_for_w.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
integrate_for_w.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
integrate_for_w.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
integrate_for_w.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
integrate_for_w.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
integrate_for_w.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
integrate_for_w.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
integrate_for_w.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
integrate_for_w.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
integrate_for_w.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
integrate_for_w.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
lef_zero.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
lef_zero.f: EEPARAMS.h SIZE.h
load_fields_driver.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
load_fields_driver.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
load_fields_driver.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
load_fields_driver.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
load_fields_driver.f: FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
load_fields_driver.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
load_fields_driver.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
load_fields_driver.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
load_fields_driver.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
load_fields_driver.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
load_fields_driver.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
load_fields_driver.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
load_fields_driver.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
load_fields_driver.f: YC_MACROS.h
load_grid_spacing.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
load_grid_spacing.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
load_grid_spacing.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
load_grid_spacing.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
load_grid_spacing.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
load_grid_spacing.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
load_grid_spacing.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
load_grid_spacing.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
load_grid_spacing.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
load_grid_spacing.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
load_grid_spacing.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
load_grid_spacing.f: RECIP_HFACW_MACROS.h SET_GRID.h SIZE.h TANPHIATU_MACROS.h
load_grid_spacing.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
load_ref_files.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
load_ref_files.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
load_ref_files.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
load_ref_files.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
load_ref_files.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
load_ref_files.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
load_ref_files.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
load_ref_files.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
load_ref_files.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
load_ref_files.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
load_ref_files.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
load_ref_files.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
main.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
main.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
main.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
main.f: EEPARAMS.h EESUPPORT.h FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
main.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MAIN_PDIRECTIVES1.h
main.f: MAIN_PDIRECTIVES2.h MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
main.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
main.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
main.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
main.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
main.f: RECIP_HFACW_MACROS.h SIGREG.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
main.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h mpif.h
main_do_loop.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
main_do_loop.f: CPP_OPTIONS.h EEPARAMS.h GAD_OPTIONS.h PACKAGES_CONFIG.h
main_do_loop.f: PARAMS.h SIZE.h
master_cpu_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
master_cpu_io.f: EEPARAMS.h SIZE.h
master_cpu_thread.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
master_cpu_thread.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
mds_byteswapi4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mds_byteswapr4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mds_byteswapr8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mds_flush.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mds_flush.f: EEPARAMS.h
mds_reclen.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mds_reclen.f: EEPARAMS.h SIZE.h
mdsfindunit.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsfindunit.f: EEPARAMS.h
mdsio_buffertorl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_buffertors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_check4file.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_check4file.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_check4file.f: PARAMS.h SIZE.h
mdsio_facef_read.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_facef_read.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_facef_read.f: SIZE.h
mdsio_gl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_gl.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_gl.f: PARAMS.h SIZE.h
mdsio_gl_slice.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_gl_slice.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_gl_slice.f: PARAMS.h SIZE.h
mdsio_pass_r4torl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_pass_r4torl.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_pass_r4torl.f: PACKAGES_CONFIG.h SIZE.h
mdsio_pass_r4tors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_pass_r4tors.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_pass_r4tors.f: PACKAGES_CONFIG.h SIZE.h
mdsio_pass_r8torl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_pass_r8torl.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_pass_r8torl.f: PACKAGES_CONFIG.h SIZE.h
mdsio_pass_r8tors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_pass_r8tors.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_pass_r8tors.f: PACKAGES_CONFIG.h SIZE.h
mdsio_rd_rec_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_rd_rec_rl.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_rd_rec_rl.f: PARAMS.h SIZE.h
mdsio_rd_rec_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_rd_rec_rs.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_rd_rec_rs.f: PARAMS.h SIZE.h
mdsio_read_field.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_read_field.f: CPP_OPTIONS.h EEBUFF_SCPU.h EEPARAMS.h MDSIO_BUFF_3D.h
mdsio_read_field.f: MDSIO_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_read_meta.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_read_meta.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_read_meta.f: PARAMS.h SIZE.h
mdsio_read_section.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mdsio_read_section.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_read_section.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_read_tape.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_read_tape.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_read_tape.f: PARAMS.h SIZE.h
mdsio_read_whalos.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_read_whalos.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_read_whalos.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_readvec_loc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_readvec_loc.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_BUFF_3D.h MDSIO_OPTIONS.h
mdsio_readvec_loc.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_rw_field.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_rw_field.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_rw_slice.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_rw_slice.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_seg4torl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_seg4torl.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_seg4tors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_seg4tors.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_seg8torl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_seg8torl.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_seg8tors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_seg8tors.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_segxtorx_2d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_segxtorx_2d.f: CPP_OPTIONS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mdsio_wr_metafiles.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mdsio_wr_metafiles.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_wr_metafiles.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_wr_rec_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_wr_rec_rl.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_wr_rec_rl.f: PARAMS.h SIZE.h
mdsio_wr_rec_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_wr_rec_rs.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_wr_rec_rs.f: PARAMS.h SIZE.h
mdsio_write_field.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_write_field.f: CPP_OPTIONS.h EEBUFF_SCPU.h EEPARAMS.h MDSIO_BUFF_3D.h
mdsio_write_field.f: MDSIO_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_write_meta.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_write_meta.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_write_meta.f: SIZE.h
mdsio_write_section.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mdsio_write_section.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mdsio_write_section.f: MDSIO_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_write_tape.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_write_tape.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h PACKAGES_CONFIG.h
mdsio_write_tape.f: PARAMS.h SIZE.h
mdsio_write_whalos.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mdsio_write_whalos.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MDSIO_OPTIONS.h
mdsio_write_whalos.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_writelocal.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mdsio_writelocal.f: CPP_OPTIONS.h EEPARAMS.h MDSIO_BUFF_3D.h MDSIO_OPTIONS.h
mdsio_writelocal.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mdsio_writevec_loc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mdsio_writevec_loc.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MDSIO_BUFF_3D.h
mdsio_writevec_loc.f: MDSIO_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
memsync.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_3d_strain.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_calc_3d_strain.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_calc_3d_strain.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_calc_3d_strain.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
mom_calc_3d_strain.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_calc_3d_strain.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_calc_3d_strain.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_calc_3d_strain.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_calc_3d_strain.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_calc_3d_strain.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_calc_3d_strain.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_calc_3d_strain.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_calc_3d_strain.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_calc_3d_strain.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_calc_absvort3.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_absvort3.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_absvort3.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_absvort3.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_calc_absvort3.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_calc_absvort3.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_calc_absvort3.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_calc_absvort3.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_calc_absvort3.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_calc_absvort3.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_calc_absvort3.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_absvort3.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_calc_absvort3.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_absvort3.f: YC_MACROS.h
mom_calc_hdiv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_hdiv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_hdiv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_hdiv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_calc_hdiv.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_calc_hdiv.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_calc_hdiv.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_calc_hdiv.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_calc_hdiv.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_calc_hdiv.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_hdiv.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_calc_hdiv.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_calc_hfacz.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_hfacz.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_hfacz.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_hfacz.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_calc_hfacz.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_calc_hfacz.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_calc_hfacz.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_calc_hfacz.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_calc_hfacz.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_calc_hfacz.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_hfacz.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_calc_hfacz.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_hfacz.f: YC_MACROS.h
mom_calc_ke.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_ke.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_ke.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_ke.f: DYU_MACROS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_calc_ke.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_calc_ke.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h RAS_MACROS.h
mom_calc_ke.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_calc_ke.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_calc_ke.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_calc_ke.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
mom_calc_ke.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_ke.f: YC_MACROS.h
mom_calc_relvort3.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_relvort3.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_relvort3.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_relvort3.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_calc_relvort3.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_calc_relvort3.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_calc_relvort3.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_calc_relvort3.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_calc_relvort3.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_calc_relvort3.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_calc_relvort3.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_relvort3.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_calc_relvort3.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_relvort3.f: YC_MACROS.h
mom_calc_rtrans.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_rtrans.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_rtrans.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_rtrans.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_calc_rtrans.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_calc_rtrans.f: MASKS_MACROS.h MASKW_MACROS.h MOM_FLUXFORM.h
mom_calc_rtrans.f: MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_calc_rtrans.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_calc_rtrans.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_calc_rtrans.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_calc_rtrans.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_rtrans.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_calc_rtrans.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_rtrans.f: YC_MACROS.h
mom_calc_smag_3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_smag_3d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_smag_3d.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_smag_3d.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_calc_smag_3d.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_calc_smag_3d.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
mom_calc_smag_3d.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_calc_smag_3d.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_calc_smag_3d.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_calc_smag_3d.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_calc_smag_3d.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_calc_smag_3d.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_calc_smag_3d.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_calc_strain.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_strain.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_strain.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_strain.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_calc_strain.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_calc_strain.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_calc_strain.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_calc_strain.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_calc_strain.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_calc_strain.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_calc_strain.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_calc_strain.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_calc_strain.f: YC_MACROS.h
mom_calc_tension.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_tension.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_tension.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_tension.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_calc_tension.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_calc_tension.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
mom_calc_tension.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_calc_tension.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_calc_tension.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_calc_tension.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_calc_tension.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_calc_tension.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_calc_tension.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_calc_visc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_calc_visc.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_calc_visc.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_calc_visc.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_calc_visc.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_calc_visc.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_VISC.h
mom_calc_visc.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_calc_visc.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_calc_visc.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_calc_visc.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_calc_visc.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_calc_visc.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_calc_visc.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_diagnostics_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_diagnostics_init.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mom_diagnostics_init.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mom_fluxform.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_fluxform.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_fluxform.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
mom_fluxform.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
mom_fluxform.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_fluxform.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
mom_fluxform.f: MOM_FLUXFORM_OPTIONS.h MOM_VISC.h PACKAGES_CONFIG.h PARAMS.h
mom_fluxform.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_fluxform.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_fluxform.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_fluxform.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_fluxform.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
mom_fluxform.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_hdissip.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_hdissip.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_hdissip.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_hdissip.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_hdissip.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_hdissip.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_hdissip.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_hdissip.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_hdissip.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_hdissip.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_hdissip.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_hdissip.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_init_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_init_fixed.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_init_fixed.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_init_fixed.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_init_fixed.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_init_fixed.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_VISC.h
mom_init_fixed.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_init_fixed.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_init_fixed.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_init_fixed.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_init_fixed.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_init_fixed.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
mom_init_fixed.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_quasihydrostatic.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_quasihydrostatic.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_quasihydrostatic.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_quasihydrostatic.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mom_quasihydrostatic.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_quasihydrostatic.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_quasihydrostatic.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_quasihydrostatic.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_quasihydrostatic.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_quasihydrostatic.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_quasihydrostatic.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_quasihydrostatic.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_quasihydrostatic.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_quasihydrostatic.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_quasihydrostatic.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_adv_uu.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_adv_uu.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_adv_uu.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_adv_uu.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_adv_uu.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_adv_uu.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_u_adv_uu.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_adv_uu.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_adv_uu.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_adv_uu.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_adv_uu.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_adv_uu.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_adv_vu.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_adv_vu.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_adv_vu.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_adv_vu.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_adv_vu.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_adv_vu.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_u_adv_vu.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_adv_vu.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_adv_vu.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_adv_vu.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_adv_vu.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_adv_vu.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_adv_wu.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_adv_wu.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_adv_wu.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_adv_wu.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_adv_wu.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_adv_wu.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_u_adv_wu.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_adv_wu.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_adv_wu.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_adv_wu.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_adv_wu.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_adv_wu.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_botdrag_coeff.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_u_botdrag_coeff.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_u_botdrag_coeff.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_u_botdrag_coeff.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_u_botdrag_coeff.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_u_botdrag_coeff.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_u_botdrag_coeff.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_u_botdrag_coeff.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_botdrag_coeff.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_botdrag_coeff.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_u_botdrag_coeff.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_botdrag_coeff.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_botdrag_coeff.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_u_botdrag_coeff.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_u_botdrag_coeff.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_coriolis.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_coriolis.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_coriolis.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_coriolis.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_coriolis.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_coriolis.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_u_coriolis.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_coriolis.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_coriolis.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_coriolis.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_coriolis.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_u_coriolis.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_coriolis.f: YC_MACROS.h
mom_u_coriolis_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_coriolis_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_coriolis_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_coriolis_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_u_coriolis_nh.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_u_coriolis_nh.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_u_coriolis_nh.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_u_coriolis_nh.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_u_coriolis_nh.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_u_coriolis_nh.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_u_coriolis_nh.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_coriolis_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_coriolis_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_coriolis_nh.f: YC_MACROS.h
mom_u_del2u.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_del2u.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_del2u.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_del2u.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_del2u.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_del2u.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_u_del2u.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_u_del2u.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_del2u.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_u_del2u.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_del2u.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
mom_u_del2u.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
mom_u_del2u.f: XC_MACROS.h YC_MACROS.h
mom_u_implicit_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_implicit_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_implicit_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_implicit_r.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
mom_u_implicit_r.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_u_implicit_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_u_implicit_r.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_u_implicit_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_u_implicit_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_u_implicit_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_u_implicit_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_implicit_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_implicit_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_implicit_r.f: YC_MACROS.h
mom_u_metric_cylinder.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_u_metric_cylinder.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_u_metric_cylinder.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_u_metric_cylinder.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mom_u_metric_cylinder.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_u_metric_cylinder.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_metric_cylinder.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h
mom_u_metric_cylinder.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_u_metric_cylinder.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_metric_cylinder.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_u_metric_cylinder.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_metric_cylinder.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_metric_cylinder.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_u_metric_cylinder.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_u_metric_cylinder.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_metric_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_metric_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_metric_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_metric_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_metric_nh.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_metric_nh.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_u_metric_nh.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_metric_nh.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_metric_nh.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_metric_nh.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_metric_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_metric_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_metric_nh.f: YC_MACROS.h
mom_u_metric_sphere.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_u_metric_sphere.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_u_metric_sphere.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_u_metric_sphere.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_u_metric_sphere.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_u_metric_sphere.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_u_metric_sphere.f: MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_u_metric_sphere.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_metric_sphere.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_metric_sphere.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_u_metric_sphere.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_metric_sphere.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_metric_sphere.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_u_metric_sphere.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_u_metric_sphere.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_rviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_rviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_rviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_rviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_rviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_rviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_u_rviscflux.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_rviscflux.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_rviscflux.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_rviscflux.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_rviscflux.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_u_rviscflux.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_rviscflux.f: YC_MACROS.h
mom_u_sidedrag.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_sidedrag.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_sidedrag.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_sidedrag.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_sidedrag.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_sidedrag.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_u_sidedrag.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_u_sidedrag.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_u_sidedrag.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_u_sidedrag.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_u_sidedrag.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_u_sidedrag.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_u_sidedrag.f: YC_MACROS.h
mom_u_xviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_xviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_xviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_xviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_xviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_xviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_u_xviscflux.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_u_xviscflux.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_xviscflux.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_u_xviscflux.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_xviscflux.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_u_xviscflux.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_u_xviscflux.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_u_yviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_u_yviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_u_yviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_u_yviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_u_yviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_u_yviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_u_yviscflux.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_u_yviscflux.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_u_yviscflux.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_u_yviscflux.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_u_yviscflux.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_u_yviscflux.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_u_yviscflux.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_uv_boundary.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_uv_boundary.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_uv_boundary.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_uv_boundary.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_uv_boundary.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_uv_boundary.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_uv_boundary.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_uv_boundary.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_uv_boundary.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_uv_boundary.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_uv_boundary.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_uv_boundary.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_uv_boundary.f: YC_MACROS.h
mom_uv_smag_3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_uv_smag_3d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_uv_smag_3d.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_uv_smag_3d.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_uv_smag_3d.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_uv_smag_3d.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_uv_smag_3d.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_uv_smag_3d.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_uv_smag_3d.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_uv_smag_3d.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_uv_smag_3d.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_uv_smag_3d.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_uv_smag_3d.f: YC_MACROS.h
mom_v_adv_uv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_adv_uv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_adv_uv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_adv_uv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_adv_uv.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_adv_uv.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_v_adv_uv.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_adv_uv.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_adv_uv.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_adv_uv.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_adv_uv.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_adv_uv.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_adv_vv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_adv_vv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_adv_vv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_adv_vv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_adv_vv.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_adv_vv.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_v_adv_vv.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_adv_vv.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_adv_vv.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_adv_vv.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_adv_vv.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_adv_vv.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_adv_wv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_adv_wv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_adv_wv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_adv_wv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_adv_wv.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_adv_wv.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_v_adv_wv.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_adv_wv.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_adv_wv.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_adv_wv.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_adv_wv.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_adv_wv.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_botdrag_coeff.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_v_botdrag_coeff.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_v_botdrag_coeff.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_v_botdrag_coeff.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_v_botdrag_coeff.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_v_botdrag_coeff.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_v_botdrag_coeff.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_v_botdrag_coeff.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_botdrag_coeff.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_botdrag_coeff.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_v_botdrag_coeff.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_botdrag_coeff.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_botdrag_coeff.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_v_botdrag_coeff.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_v_botdrag_coeff.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_coriolis.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_coriolis.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_coriolis.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_coriolis.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_coriolis.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_coriolis.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h
mom_v_coriolis.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_coriolis.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_coriolis.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_coriolis.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_coriolis.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_v_coriolis.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_coriolis.f: YC_MACROS.h
mom_v_coriolis_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_coriolis_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_coriolis_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_coriolis_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_v_coriolis_nh.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_v_coriolis_nh.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_v_coriolis_nh.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_v_coriolis_nh.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_v_coriolis_nh.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_v_coriolis_nh.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_v_coriolis_nh.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_coriolis_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_coriolis_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_coriolis_nh.f: YC_MACROS.h
mom_v_del2v.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_del2v.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_del2v.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_del2v.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_del2v.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_del2v.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_v_del2v.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_v_del2v.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_del2v.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_v_del2v.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_del2v.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
mom_v_del2v.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
mom_v_del2v.f: XC_MACROS.h YC_MACROS.h
mom_v_implicit_r.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_implicit_r.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_implicit_r.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_implicit_r.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
mom_v_implicit_r.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_v_implicit_r.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_v_implicit_r.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_v_implicit_r.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_v_implicit_r.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_v_implicit_r.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_v_implicit_r.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_implicit_r.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_implicit_r.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_implicit_r.f: YC_MACROS.h
mom_v_metric_cylinder.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_v_metric_cylinder.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_v_metric_cylinder.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_v_metric_cylinder.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mom_v_metric_cylinder.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_v_metric_cylinder.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_metric_cylinder.f: MASKW_MACROS.h MOM_FLUXFORM_OPTIONS.h
mom_v_metric_cylinder.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_v_metric_cylinder.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_metric_cylinder.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_v_metric_cylinder.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_metric_cylinder.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_metric_cylinder.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_v_metric_cylinder.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_v_metric_cylinder.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_metric_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_metric_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_metric_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_metric_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_metric_nh.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_metric_nh.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_v_metric_nh.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_metric_nh.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_metric_nh.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_metric_nh.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_metric_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_metric_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_metric_nh.f: YC_MACROS.h
mom_v_metric_sphere.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_v_metric_sphere.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_v_metric_sphere.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_v_metric_sphere.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_v_metric_sphere.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_v_metric_sphere.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_v_metric_sphere.f: MOM_FLUXFORM_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_v_metric_sphere.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_metric_sphere.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_metric_sphere.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_v_metric_sphere.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_metric_sphere.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_metric_sphere.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_v_metric_sphere.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_v_metric_sphere.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_rviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_rviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_rviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_rviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_rviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_rviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_v_rviscflux.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_rviscflux.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_rviscflux.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_rviscflux.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_rviscflux.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_v_rviscflux.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_rviscflux.f: YC_MACROS.h
mom_v_sidedrag.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_sidedrag.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_sidedrag.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_sidedrag.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_sidedrag.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_sidedrag.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_v_sidedrag.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_v_sidedrag.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_v_sidedrag.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_v_sidedrag.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_v_sidedrag.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_v_sidedrag.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_v_sidedrag.f: YC_MACROS.h
mom_v_xviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_xviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_xviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_xviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_xviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_xviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_v_xviscflux.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_v_xviscflux.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_xviscflux.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_v_xviscflux.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_xviscflux.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_v_xviscflux.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_v_xviscflux.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_v_yviscflux.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_v_yviscflux.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_v_yviscflux.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_v_yviscflux.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_v_yviscflux.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_v_yviscflux.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h MOM_FLUXFORM_OPTIONS.h
mom_v_yviscflux.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_v_yviscflux.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_v_yviscflux.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mom_v_yviscflux.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_v_yviscflux.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_v_yviscflux.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_v_yviscflux.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_vecinv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vecinv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vecinv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
mom_vecinv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
mom_vecinv.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_vecinv.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
mom_vecinv.f: MOM_VECINV_OPTIONS.h MOM_VISC.h PACKAGES_CONFIG.h PARAMS.h
mom_vecinv.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vecinv.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vecinv.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vecinv.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_vecinv.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
mom_vecinv.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_vi_coriolis.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_coriolis.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_coriolis.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_coriolis.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_vi_coriolis.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_vi_coriolis.f: MASKW_MACROS.h MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h
mom_vi_coriolis.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_coriolis.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_vi_coriolis.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_coriolis.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_coriolis.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_coriolis.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_coriolis.f: YC_MACROS.h
mom_vi_del2uv.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_del2uv.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_del2uv.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_del2uv.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_vi_del2uv.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_vi_del2uv.f: MASKW_MACROS.h MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h
mom_vi_del2uv.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vi_del2uv.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_del2uv.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vi_del2uv.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_vi_del2uv.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_vi_del2uv.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_vi_hdissip.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_hdissip.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_hdissip.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_hdissip.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_vi_hdissip.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_vi_hdissip.f: MASKW_MACROS.h MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h
mom_vi_hdissip.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_hdissip.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_vi_hdissip.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_hdissip.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_hdissip.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_hdissip.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_hdissip.f: YC_MACROS.h
mom_vi_u_coriolis.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_u_coriolis.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_u_coriolis.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_u_coriolis.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_vi_u_coriolis.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_vi_u_coriolis.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_vi_u_coriolis.f: MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h RAS_MACROS.h
mom_vi_u_coriolis.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vi_u_coriolis.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_u_coriolis.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vi_u_coriolis.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_u_coriolis.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_u_coriolis.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_u_coriolis.f: YC_MACROS.h
mom_vi_u_coriolis_c4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_vi_u_coriolis_c4.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_vi_u_coriolis_c4.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_vi_u_coriolis_c4.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mom_vi_u_coriolis_c4.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_vi_u_coriolis_c4.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_vi_u_coriolis_c4.f: MASKW_MACROS.h MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h
mom_vi_u_coriolis_c4.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_vi_u_coriolis_c4.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_vi_u_coriolis_c4.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_vi_u_grad_ke.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_u_grad_ke.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_u_grad_ke.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_u_grad_ke.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_vi_u_grad_ke.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_vi_u_grad_ke.f: MASKS_MACROS.h MASKW_MACROS.h MOM_VECINV_OPTIONS.h
mom_vi_u_grad_ke.f: PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_u_grad_ke.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_vi_u_grad_ke.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_u_grad_ke.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_u_grad_ke.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_u_grad_ke.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_u_grad_ke.f: YC_MACROS.h
mom_vi_u_vertshear.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_vi_u_vertshear.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_vi_u_vertshear.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_vi_u_vertshear.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_vi_u_vertshear.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_vi_u_vertshear.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_vi_u_vertshear.f: MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_vi_u_vertshear.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vi_u_vertshear.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_u_vertshear.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vi_u_vertshear.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_u_vertshear.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_u_vertshear.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_u_vertshear.f: YC_MACROS.h
mom_vi_v_coriolis.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_v_coriolis.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_v_coriolis.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_v_coriolis.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_vi_v_coriolis.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_vi_v_coriolis.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_vi_v_coriolis.f: MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h RAS_MACROS.h
mom_vi_v_coriolis.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vi_v_coriolis.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_v_coriolis.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vi_v_coriolis.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_v_coriolis.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_v_coriolis.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_v_coriolis.f: YC_MACROS.h
mom_vi_v_coriolis_c4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_vi_v_coriolis_c4.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_vi_v_coriolis_c4.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_vi_v_coriolis_c4.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mom_vi_v_coriolis_c4.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mom_vi_v_coriolis_c4.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_vi_v_coriolis_c4.f: MASKW_MACROS.h MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h
mom_vi_v_coriolis_c4.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_vi_v_coriolis_c4.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_vi_v_coriolis_c4.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_vi_v_grad_ke.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_vi_v_grad_ke.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_vi_v_grad_ke.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_vi_v_grad_ke.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_vi_v_grad_ke.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_vi_v_grad_ke.f: MASKS_MACROS.h MASKW_MACROS.h MOM_VECINV_OPTIONS.h
mom_vi_v_grad_ke.f: PACKAGES_CONFIG.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_vi_v_grad_ke.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_vi_v_grad_ke.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_vi_v_grad_ke.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_v_grad_ke.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_v_grad_ke.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_v_grad_ke.f: YC_MACROS.h
mom_vi_v_vertshear.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_vi_v_vertshear.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_vi_v_vertshear.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_vi_v_vertshear.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_vi_v_vertshear.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_vi_v_vertshear.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_vi_v_vertshear.f: MOM_VECINV_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_vi_v_vertshear.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_vi_v_vertshear.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_vi_v_vertshear.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_vi_v_vertshear.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_vi_v_vertshear.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_vi_v_vertshear.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_vi_v_vertshear.f: YC_MACROS.h
mom_visc_qgl_limit.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_visc_qgl_limit.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
mom_visc_qgl_limit.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
mom_visc_qgl_limit.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h
mom_visc_qgl_limit.f: GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_visc_qgl_limit.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_visc_qgl_limit.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_visc_qgl_limit.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_visc_qgl_limit.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_visc_qgl_limit.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_visc_qgl_limit.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_visc_qgl_limit.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_visc_qgl_limit.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_visc_qgl_limit.f: YC_MACROS.h
mom_visc_qgl_stretch.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mom_visc_qgl_stretch.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mom_visc_qgl_stretch.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mom_visc_qgl_stretch.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
mom_visc_qgl_stretch.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_visc_qgl_stretch.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
mom_visc_qgl_stretch.f: MASKS_MACROS.h MASKW_MACROS.h MOM_COMMON_OPTIONS.h
mom_visc_qgl_stretch.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mom_visc_qgl_stretch.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mom_visc_qgl_stretch.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_visc_qgl_stretch.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_visc_qgl_stretch.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mom_visc_qgl_stretch.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mom_visc_qgl_stretch.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mom_visc_qgl_stretch.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mom_w_coriolis_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_w_coriolis_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_w_coriolis_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_w_coriolis_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mom_w_coriolis_nh.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mom_w_coriolis_nh.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mom_w_coriolis_nh.f: MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mom_w_coriolis_nh.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mom_w_coriolis_nh.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mom_w_coriolis_nh.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mom_w_coriolis_nh.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_w_coriolis_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_w_coriolis_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_w_coriolis_nh.f: YC_MACROS.h
mom_w_metric_nh.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_w_metric_nh.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_w_metric_nh.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_w_metric_nh.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_w_metric_nh.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_w_metric_nh.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_w_metric_nh.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_w_metric_nh.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_w_metric_nh.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_w_metric_nh.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_w_metric_nh.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_w_metric_nh.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_w_metric_nh.f: YC_MACROS.h
mom_w_sidedrag.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_w_sidedrag.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_w_sidedrag.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_w_sidedrag.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_w_sidedrag.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_w_sidedrag.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_w_sidedrag.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_w_sidedrag.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_w_sidedrag.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_w_sidedrag.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_w_sidedrag.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mom_w_sidedrag.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mom_w_sidedrag.f: YC_MACROS.h
mom_w_smag_3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mom_w_smag_3d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mom_w_smag_3d.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mom_w_smag_3d.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mom_w_smag_3d.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mom_w_smag_3d.f: MASKW_MACROS.h MOM_COMMON_OPTIONS.h PACKAGES_CONFIG.h
mom_w_smag_3d.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mom_w_smag_3d.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mom_w_smag_3d.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mom_w_smag_3d.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mom_w_smag_3d.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mom_w_smag_3d.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
momentum_correction_step.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
momentum_correction_step.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DYNVARS.h EEPARAMS.h
momentum_correction_step.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mon_advcfl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_advcfl.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_advcfl.f: PACKAGES_CONFIG.h SIZE.h
mon_advcflw.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_advcflw.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_advcflw.f: PACKAGES_CONFIG.h SIZE.h
mon_advcflw2.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_advcflw2.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_advcflw2.f: PACKAGES_CONFIG.h SIZE.h
mon_calc_advcfl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_calc_advcfl.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mon_calc_advcfl.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mon_calc_advcfl.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mon_calc_advcfl.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mon_calc_advcfl.f: MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h
mon_calc_advcfl.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mon_calc_advcfl.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mon_calc_advcfl.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mon_calc_advcfl.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mon_calc_advcfl.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mon_calc_advcfl.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mon_calc_advcfl.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mon_calc_stats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_calc_stats_rl.f: CPP_OPTIONS.h EEPARAMS.h MONITOR_OPTIONS.h
mon_calc_stats_rl.f: PACKAGES_CONFIG.h SIZE.h
mon_calc_stats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_calc_stats_rs.f: CPP_OPTIONS.h EEPARAMS.h MONITOR_OPTIONS.h
mon_calc_stats_rs.f: PACKAGES_CONFIG.h SIZE.h
mon_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_init.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_init.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mon_ke.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_ke.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
mon_ke.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
mon_ke.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mon_ke.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mon_ke.f: MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
mon_ke.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mon_ke.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
mon_ke.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mon_ke.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
mon_ke.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
mon_ke.f: YC_MACROS.h
mon_out.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_out.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_out.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mon_printstats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_printstats_rl.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_printstats_rl.f: PACKAGES_CONFIG.h SIZE.h
mon_printstats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_printstats_rs.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_printstats_rs.f: PACKAGES_CONFIG.h SIZE.h
mon_set_iounit.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_set_iounit.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_set_iounit.f: PACKAGES_CONFIG.h SIZE.h
mon_set_pref.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_set_pref.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_set_pref.f: PACKAGES_CONFIG.h SIZE.h
mon_solution.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_solution.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mon_solution.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
mon_solution.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mon_solution.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mon_solution.f: MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
mon_solution.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mon_solution.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mon_solution.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mon_solution.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mon_solution.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
mon_solution.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mon_stats_latbnd_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mon_stats_latbnd_rl.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mon_stats_latbnd_rl.f: MONITOR_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mon_stats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_stats_rl.f: CPP_OPTIONS.h EEPARAMS.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
mon_stats_rl.f: SIZE.h
mon_stats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_stats_rs.f: CPP_OPTIONS.h EEPARAMS.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
mon_stats_rs.f: SIZE.h
mon_surfcor.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_surfcor.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
mon_surfcor.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
mon_surfcor.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
mon_surfcor.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mon_surfcor.f: MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
mon_surfcor.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mon_surfcor.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
mon_surfcor.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mon_surfcor.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
mon_surfcor.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
mon_surfcor.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mon_vort3.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_vort3.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
mon_vort3.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
mon_vort3.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mon_vort3.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mon_vort3.f: MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mon_vort3.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
mon_vort3.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mon_vort3.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
mon_vort3.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mon_vort3.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
mon_vort3.f: XC_MACROS.h YC_MACROS.h
mon_writestats_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_writestats_rl.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_writestats_rl.f: PACKAGES_CONFIG.h SIZE.h
mon_writestats_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mon_writestats_rs.f: CPP_OPTIONS.h EEPARAMS.h MONITOR.h MONITOR_OPTIONS.h
mon_writestats_rs.f: PACKAGES_CONFIG.h SIZE.h
monitor.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
monitor.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
monitor.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
monitor.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
monitor.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
monitor.f: MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
monitor.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
monitor.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
monitor.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
monitor.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
monitor.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
monitor.f: XC_MACROS.h YC_MACROS.h
monitor_ad.f: /usr/include/stdc-predef.h AD_CONFIG.h CPP_EEMACROS.h
monitor_ad.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
monitor_ad.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
monitor_ad.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
monitor_ad.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
monitor_ad.f: MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h PACKAGES_CONFIG.h
monitor_ad.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
monitor_ad.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
monitor_ad.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
monitor_ad.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
monitor_ad.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
monitor_ad.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
monitor_g.f: /usr/include/stdc-predef.h AD_CONFIG.h CPP_EEMACROS.h
monitor_g.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
monitor_g.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
monitor_g.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h
monitor_g.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
monitor_g.f: MASKS_MACROS.h MASKW_MACROS.h MONITOR.h MONITOR_OPTIONS.h
monitor_g.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
monitor_g.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
monitor_g.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
monitor_g.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
monitor_g.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
monitor_g.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mypackage_calc_rhs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_calc_rhs.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mypackage_calc_rhs.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h SIZE.h
mypackage_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
mypackage_check.f: CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h MYPACKAGE_OPTIONS.h
mypackage_check.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_diagnostics_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_diagnostics_init.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mypackage_diagnostics_init.f: MYPACKAGE.h MYPACKAGE_OPTIONS.h
mypackage_diagnostics_init.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_diagnostics_state.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_diagnostics_state.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
mypackage_diagnostics_state.f: MYPACKAGE.h MYPACKAGE_OPTIONS.h
mypackage_diagnostics_state.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_init_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_init_fixed.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h
mypackage_init_fixed.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_init_varia.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_init_varia.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mypackage_init_varia.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
mypackage_init_varia.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
mypackage_init_varia.f: FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
mypackage_init_varia.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
mypackage_init_varia.f: MASKW_MACROS.h MYPACKAGE.h MYPACKAGE_OPTIONS.h
mypackage_init_varia.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
mypackage_init_varia.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mypackage_init_varia.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mypackage_init_varia.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mypackage_init_varia.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mypackage_init_varia.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mypackage_init_varia.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mypackage_init_varia.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mypackage_mnc_init.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_mnc_init.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h
mypackage_mnc_init.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_read_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_read_pickup.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h
mypackage_read_pickup.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_readparms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_readparms.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h
mypackage_readparms.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
mypackage_tendency_apply.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_tendency_apply.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
mypackage_tendency_apply.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
mypackage_tendency_apply.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
mypackage_tendency_apply.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
mypackage_tendency_apply.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
mypackage_tendency_apply.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
mypackage_tendency_apply.f: MYPACKAGE.h MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h
mypackage_tendency_apply.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
mypackage_tendency_apply.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
mypackage_tendency_apply.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
mypackage_tendency_apply.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
mypackage_tendency_apply.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
mypackage_tendency_apply.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
mypackage_tendency_apply.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
mypackage_tendency_apply.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
mypackage_write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
mypackage_write_pickup.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h MYPACKAGE.h
mypackage_write_pickup.f: MYPACKAGE_OPTIONS.h PACKAGES_CONFIG.h PARAMS.h
mypackage_write_pickup.f: SIZE.h
nml_change_syntax.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
nml_change_syntax.f: EEPARAMS.h
nml_set_terminator.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
nml_set_terminator.f: CPP_EEOPTIONS.h EEPARAMS.h
open_copy_data_file.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
open_copy_data_file.f: CPP_EEOPTIONS.h EEPARAMS.h SIZE.h
packages_boot.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
packages_boot.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
packages_check.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
packages_check.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
packages_error_msg.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_error_msg.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h SIZE.h
packages_init_fixed.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_init_fixed.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
packages_init_fixed.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
packages_init_variables.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_init_variables.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
packages_init_variables.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
packages_print_msg.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_print_msg.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h SIZE.h
packages_readparms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_readparms.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
packages_readparms.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
packages_unused_msg.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_unused_msg.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h SIZE.h
packages_write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
packages_write_pickup.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
packages_write_pickup.f: PACKAGES_CONFIG.h PARAMS.h RESTART.h SIZE.h
plot_field.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
plot_field.f: CPP_OPTIONS.h EEPARAMS.h PARAMS.h SIZE.h
port_rand.f: /usr/include/stdc-predef.h
post_cg3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
post_cg3d.f: CPP_OPTIONS.h EEPARAMS.h SIZE.h
pre_cg3d.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
pre_cg3d.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
pre_cg3d.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
pre_cg3d.f: EEPARAMS.h FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
pre_cg3d.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
pre_cg3d.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
pre_cg3d.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
pre_cg3d.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
pre_cg3d.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
pre_cg3d.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
pre_cg3d.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
pre_cg3d.f: YC_MACROS.h
pressure_for_eos.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
pressure_for_eos.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
pressure_for_eos.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
pressure_for_eos.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
pressure_for_eos.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
pressure_for_eos.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
pressure_for_eos.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
pressure_for_eos.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
pressure_for_eos.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
pressure_for_eos.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
pressure_for_eos.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
pressure_for_eos.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
pressure_for_eos.f: YC_MACROS.h
print.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h EEPARAMS.h
print.f: EESUPPORT.h SIZE.h mpif.h
read_fld_xy_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_fld_xy_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_fld_xy_rl.f: RW_OPTIONS.h SIZE.h
read_fld_xy_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_fld_xy_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_fld_xy_rs.f: RW_OPTIONS.h SIZE.h
read_fld_xyz_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_fld_xyz_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_fld_xyz_rl.f: RW_OPTIONS.h SIZE.h
read_fld_xyz_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_fld_xyz_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_fld_xyz_rs.f: RW_OPTIONS.h SIZE.h
read_glvec_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_glvec_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_glvec_rl.f: RW_OPTIONS.h SIZE.h
read_glvec_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_glvec_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
read_glvec_rs.f: RW_OPTIONS.h SIZE.h
read_mflds.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_mflds.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h RW_MFLDS.h
read_mflds.f: RW_OPTIONS.h SIZE.h
read_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_pickup.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h FFIELDS.h GAD.h NH_VARS.h
read_pickup.f: PACKAGES_CONFIG.h PARAMS.h RESTART.h SIZE.h SURFACE.h
read_rec.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
read_rec.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h RW_OPTIONS.h
read_rec.f: SIZE.h
remove_mean.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
remove_mean.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
remove_mean.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
remove_mean.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
remove_mean.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
remove_mean.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
remove_mean.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
remove_mean.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
remove_mean.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
remove_mean.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
remove_mean.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
reset_halo.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
reset_halo.f: EEPARAMS.h SIZE.h
reset_nlfs_vars.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
reset_nlfs_vars.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
reset_nlfs_vars.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
reset_nlfs_vars.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
reset_nlfs_vars.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
reset_nlfs_vars.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
reset_nlfs_vars.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
reset_nlfs_vars.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
reset_nlfs_vars.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
reset_nlfs_vars.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
reset_nlfs_vars.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
reset_nlfs_vars.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
rotate_spherical_polar_grid.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
rotate_spherical_polar_grid.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
rotate_spherical_polar_grid.f: PARAMS.h SIZE.h
rotate_uv2en.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
rotate_uv2en.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
rotate_uv2en.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
rotate_uv2en.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
rotate_uv2en.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
rotate_uv2en.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
rotate_uv2en.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
rotate_uv2en.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
rotate_uv2en.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
rotate_uv2en.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
rotate_uv2en.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
rw_get_suffix.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
rw_get_suffix.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
rw_get_suffix.f: RW_OPTIONS.h SIZE.h
salt_fill.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
salt_fill.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
salt_fill.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
salt_fill.f: FCORI_MACROS.h GAD_OPTIONS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
salt_fill.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
salt_fill.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
salt_fill.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
salt_fill.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
salt_fill.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
salt_fill.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
salt_fill.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
salt_integrate.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
salt_integrate.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
salt_integrate.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
salt_integrate.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h
salt_integrate.f: GAD_OPTIONS.h GAD_SOM_VARS.h GRID.h GRID_MACROS.h
salt_integrate.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
salt_integrate.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
salt_integrate.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
salt_integrate.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
salt_integrate.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
salt_integrate.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
salt_integrate.f: RECIP_HFACW_MACROS.h RESTART.h SIZE.h TANPHIATU_MACROS.h
salt_integrate.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
scatter_2d_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_2d_r4.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
scatter_2d_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_2d_r8.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
scatter_2d_wh_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_2d_wh_r4.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
scatter_2d_wh_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_2d_wh_r8.f: EEPARAMS.h EESUPPORT.h PACKAGES_CONFIG.h SIZE.h mpif.h
scatter_vec_r4.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_vec_r4.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
scatter_vec_r8.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_vec_r8.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
scatter_xz.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_xz.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
scatter_yz.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
scatter_yz.f: EEPARAMS.h EESUPPORT.h SIZE.h mpif.h
seawater.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
seawater.f: CPP_OPTIONS.h EOS.h PACKAGES_CONFIG.h SIZE.h
set_defaults.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
set_defaults.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SET_GRID.h
set_defaults.f: SIZE.h
set_grid_factors.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
set_grid_factors.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
set_grid_factors.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
set_grid_factors.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
set_grid_factors.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
set_grid_factors.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
set_grid_factors.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
set_grid_factors.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
set_grid_factors.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
set_grid_factors.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
set_grid_factors.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
set_grid_factors.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
set_grid_factors.f: YC_MACROS.h
set_parms.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
set_parms.f: CPP_OPTIONS.h EEPARAMS.h EOS.h MOM_COMMON_OPTIONS.h MOM_VISC.h
set_parms.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
set_ref_state.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
set_ref_state.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
set_ref_state.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
set_ref_state.f: DYU_MACROS.h EEPARAMS.h EOS.h FCORI_MACROS.h GRID.h
set_ref_state.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
set_ref_state.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
set_ref_state.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
set_ref_state.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
set_ref_state.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
set_ref_state.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
set_ref_state.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
set_ref_state.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
set_write_global_fld.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
set_write_global_fld.f: CPP_EEOPTIONS.h
solve_for_pressure.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
solve_for_pressure.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
solve_for_pressure.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
solve_for_pressure.f: DYG_MACROS.h DYNVARS.h DYU_MACROS.h EEPARAMS.h
solve_for_pressure.f: FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
solve_for_pressure.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
solve_for_pressure.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
solve_for_pressure.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
solve_for_pressure.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
solve_for_pressure.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
solve_for_pressure.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
solve_for_pressure.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
solve_for_pressure.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
solve_for_pressure.f: XC_MACROS.h YC_MACROS.h
solve_pentadiagonal.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
solve_pentadiagonal.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h
solve_pentadiagonal.f: PACKAGES_CONFIG.h SIZE.h
solve_tridiagonal.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
solve_tridiagonal.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h SIZE.h
solve_uv_tridiago.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
solve_uv_tridiago.f: CPP_OPTIONS.h EEPARAMS.h SIZE.h
state_summary.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
state_summary.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
state_summary.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
state_summary.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
state_summary.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
state_summary.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
state_summary.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
state_summary.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
state_summary.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
state_summary.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
state_summary.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
stop_if_error.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
stop_if_error.f: EEPARAMS.h SIZE.h
swfrac.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
swfrac.f: CPP_OPTIONS.h PACKAGES_CONFIG.h
taueddy_init_varia.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
taueddy_init_varia.f: CPP_EEOPTIONS.h CPP_OPTIONS.h EEPARAMS.h FFIELDS.h
taueddy_init_varia.f: PARAMS.h SIZE.h
taueddy_tendency_apply.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
taueddy_tendency_apply.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
taueddy_tendency_apply.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
taueddy_tendency_apply.f: DYF_MACROS.h DYG_MACROS.h DYU_MACROS.h EEPARAMS.h
taueddy_tendency_apply.f: FCORI_MACROS.h FFIELDS.h GRID.h GRID_MACROS.h
taueddy_tendency_apply.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
taueddy_tendency_apply.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h
taueddy_tendency_apply.f: PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
taueddy_tendency_apply.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
taueddy_tendency_apply.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
taueddy_tendency_apply.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
taueddy_tendency_apply.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
taueddy_tendency_apply.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
taueddy_tendency_apply.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
taueddy_tendency_apply.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
temp_integrate.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
temp_integrate.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
temp_integrate.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
temp_integrate.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GAD.h
temp_integrate.f: GAD_OPTIONS.h GAD_SOM_VARS.h GRID.h GRID_MACROS.h
temp_integrate.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
temp_integrate.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
temp_integrate.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
temp_integrate.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
temp_integrate.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
temp_integrate.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
temp_integrate.f: RECIP_HFACW_MACROS.h RESTART.h SIZE.h TANPHIATU_MACROS.h
temp_integrate.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
the_main_loop.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
the_main_loop.f: CPP_OPTIONS.h EEPARAMS.h GAD_OPTIONS.h PACKAGES_CONFIG.h
the_main_loop.f: PARAMS.h SIZE.h
the_model_main.f: /usr/include/stdc-predef.h AD_CONFIG.h CPP_EEMACROS.h
the_model_main.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DYNVARS.h EEPARAMS.h
the_model_main.f: PACKAGES_CONFIG.h PARAMS.h SIZE.h
thermodynamics.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
thermodynamics.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
thermodynamics.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
thermodynamics.f: DYNVARS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h FFIELDS.h
thermodynamics.f: GAD.h GRID.h GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h
thermodynamics.f: HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
thermodynamics.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
thermodynamics.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
thermodynamics.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
thermodynamics.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
thermodynamics.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
thermodynamics.f: RECIP_HFACW_MACROS.h RESTART.h SIZE.h SURFACE.h
thermodynamics.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
thermodynamics.f: YC_MACROS.h
timers.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h EEPARAMS.h
timers.f: EESUPPORT.h SIZE.h mpif.h
timestep.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
timestep.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h
timestep.f: DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
timestep.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h HFACC_MACROS.h
timestep.f: HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h MASKW_MACROS.h
timestep.f: PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
timestep.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
timestep.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
timestep.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
timestep.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h RESTART.h SIZE.h
timestep.f: SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
timestep.f: YC_MACROS.h
timestep_tracer.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
timestep_tracer.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
timestep_tracer.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
timestep_tracer.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
timestep_tracer.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
timestep_tracer.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
timestep_tracer.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
timestep_tracer.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
timestep_tracer.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
timestep_tracer.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
timestep_tracer.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
timestep_tracer.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
timestep_tracer.f: YC_MACROS.h
timestep_wvel.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
timestep_wvel.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
timestep_wvel.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
timestep_wvel.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
timestep_wvel.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
timestep_wvel.f: MASKW_MACROS.h NH_VARS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
timestep_wvel.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
timestep_wvel.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
timestep_wvel.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
timestep_wvel.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
timestep_wvel.f: RECIP_HFACW_MACROS.h SIZE.h TANPHIATU_MACROS.h
timestep_wvel.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
tracers_correction_step.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
tracers_correction_step.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DYNVARS.h EEPARAMS.h
tracers_correction_step.f: GAD.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
tracers_iigw_correction.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
tracers_iigw_correction.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h
tracers_iigw_correction.f: DXF_MACROS.h DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h
tracers_iigw_correction.f: DYF_MACROS.h DYG_MACROS.h DYNVARS.h DYU_MACROS.h
tracers_iigw_correction.f: EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
tracers_iigw_correction.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
tracers_iigw_correction.f: MASKS_MACROS.h MASKW_MACROS.h PARAMS.h RAS_MACROS.h
tracers_iigw_correction.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
tracers_iigw_correction.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
tracers_iigw_correction.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
tracers_iigw_correction.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
tracers_iigw_correction.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
tracers_iigw_correction.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h
tracers_iigw_correction.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
tracers_iigw_correction.f: YC_MACROS.h
turnoff_model_io.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
turnoff_model_io.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h SIZE.h
update_cg2d.f: /usr/include/stdc-predef.h CG2D.h CPP_EEMACROS.h
update_cg2d.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
update_cg2d.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h
update_cg2d.f: DYG_MACROS.h DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
update_cg2d.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
update_cg2d.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
update_cg2d.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
update_cg2d.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
update_cg2d.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
update_cg2d.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
update_cg2d.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
update_cg2d.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
update_etah.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_etah.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_etah.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
update_etah.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
update_etah.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
update_etah.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
update_etah.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
update_etah.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
update_etah.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
update_etah.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
update_etah.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
update_etah.f: XC_MACROS.h YC_MACROS.h
update_etaws.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_etaws.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_etaws.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
update_etaws.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
update_etaws.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
update_etaws.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
update_etaws.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
update_etaws.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
update_etaws.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
update_etaws.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
update_etaws.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
update_masks_etc.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_masks_etc.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_masks_etc.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
update_masks_etc.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h
update_masks_etc.f: GRID_MACROS.h HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h
update_masks_etc.f: MASKS_MACROS.h MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h
update_masks_etc.f: RAS_MACROS.h RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h
update_masks_etc.f: RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h
update_masks_etc.f: RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h
update_masks_etc.f: RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
update_masks_etc.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
update_masks_etc.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
update_masks_etc.f: YC_MACROS.h
update_r_star.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_r_star.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_r_star.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
update_r_star.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
update_r_star.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
update_r_star.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h RA_MACROS.h
update_r_star.f: RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h RECIP_DXG_MACROS.h
update_r_star.f: RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h RECIP_DYF_MACROS.h
update_r_star.f: RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h RECIP_HFACC_MACROS.h
update_r_star.f: RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h SIZE.h SURFACE.h
update_r_star.f: TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
update_sigma.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_sigma.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_sigma.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
update_sigma.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
update_sigma.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
update_sigma.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
update_sigma.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
update_sigma.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
update_sigma.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
update_sigma.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
update_sigma.f: SIZE.h SURFACE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h
update_sigma.f: XC_MACROS.h YC_MACROS.h
update_surf_dr.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
update_surf_dr.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
update_surf_dr.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
update_surf_dr.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
update_surf_dr.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
update_surf_dr.f: MASKW_MACROS.h PARAMS.h RAS_MACROS.h RAW_MACROS.h
update_surf_dr.f: RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
update_surf_dr.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
update_surf_dr.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
update_surf_dr.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h
update_surf_dr.f: RECIP_HFACW_MACROS.h SIZE.h SURFACE.h TANPHIATU_MACROS.h
update_surf_dr.f: TANPHIATV_MACROS.h XC_MACROS.h YC_MACROS.h
utils.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h EEPARAMS.h
utils.f: SIZE.h
write_fld_3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_3d_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_3d_rl.f: RW_OPTIONS.h SIZE.h
write_fld_3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_3d_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_3d_rs.f: RW_OPTIONS.h SIZE.h
write_fld_s3d_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_s3d_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_s3d_rl.f: RW_OPTIONS.h SIZE.h
write_fld_s3d_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_s3d_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_s3d_rs.f: RW_OPTIONS.h SIZE.h
write_fld_xy_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_xy_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_xy_rl.f: RW_OPTIONS.h SIZE.h
write_fld_xy_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_xy_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_xy_rs.f: RW_OPTIONS.h SIZE.h
write_fld_xyz_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_xyz_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_xyz_rl.f: RW_OPTIONS.h SIZE.h
write_fld_xyz_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_fld_xyz_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_fld_xyz_rs.f: RW_OPTIONS.h SIZE.h
write_fullarray_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
write_fullarray_rl.f: CPP_EEOPTIONS.h EEPARAMS.h SIZE.h
write_fullarray_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h
write_fullarray_rs.f: CPP_EEOPTIONS.h EEPARAMS.h SIZE.h
write_glvec_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_glvec_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_glvec_rl.f: RW_OPTIONS.h SIZE.h
write_glvec_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_glvec_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_glvec_rs.f: RW_OPTIONS.h SIZE.h
write_grid.f: /usr/include/stdc-predef.h AD_CONFIG.h CPP_EEMACROS.h
write_grid.f: CPP_EEOPTIONS.h CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h
write_grid.f: DXG_MACROS.h DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h
write_grid.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
write_grid.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
write_grid.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
write_grid.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
write_grid.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
write_grid.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
write_grid.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
write_grid.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
write_grid.f: YC_MACROS.h
write_local_rl.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_local_rl.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_local_rl.f: RW_OPTIONS.h SIZE.h
write_local_rs.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_local_rs.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h
write_local_rs.f: RW_OPTIONS.h SIZE.h
write_pickup.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_pickup.f: CPP_OPTIONS.h DYNVARS.h EEPARAMS.h FFIELDS.h GAD.h NH_VARS.h
write_pickup.f: PACKAGES_CONFIG.h PARAMS.h RESTART.h SIZE.h SURFACE.h
write_rec.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_rec.f: CPP_OPTIONS.h EEPARAMS.h PACKAGES_CONFIG.h PARAMS.h RW_OPTIONS.h
write_rec.f: SIZE.h
write_state.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_state.f: CPP_OPTIONS.h DXC_MACROS.h DXF_MACROS.h DXG_MACROS.h
write_state.f: DXV_MACROS.h DYC_MACROS.h DYF_MACROS.h DYG_MACROS.h DYNVARS.h
write_state.f: DYU_MACROS.h EEPARAMS.h FCORI_MACROS.h GRID.h GRID_MACROS.h
write_state.f: HFACC_MACROS.h HFACS_MACROS.h HFACW_MACROS.h MASKS_MACROS.h
write_state.f: MASKW_MACROS.h PACKAGES_CONFIG.h PARAMS.h RAS_MACROS.h
write_state.f: RAW_MACROS.h RA_MACROS.h RECIP_DXC_MACROS.h RECIP_DXF_MACROS.h
write_state.f: RECIP_DXG_MACROS.h RECIP_DXV_MACROS.h RECIP_DYC_MACROS.h
write_state.f: RECIP_DYF_MACROS.h RECIP_DYG_MACROS.h RECIP_DYU_MACROS.h
write_state.f: RECIP_HFACC_MACROS.h RECIP_HFACS_MACROS.h RECIP_HFACW_MACROS.h
write_state.f: SIZE.h TANPHIATU_MACROS.h TANPHIATV_MACROS.h XC_MACROS.h
write_state.f: YC_MACROS.h
write_utils.f: /usr/include/stdc-predef.h CPP_EEMACROS.h CPP_EEOPTIONS.h
write_utils.f: EEPARAMS.h SIZE.h

