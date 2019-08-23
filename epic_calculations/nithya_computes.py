def computations():
    fft_coeff = 5.0
    
    telescopes = ['HERA-19 (14m)', 'HERA-37 (14m)', 'HERA-331 (14m)', 'HERA-6769 (14m)', 'CHIME', 'HIRAX', 'MWA-112 (4m)', 'MWA-240 (4m)', 'MWA-496 (4m)', 'MWA-1008 (4m)', 'MWA-II-C (4m)', 'SKA1-LC (35m)', 'SKA1-LCD (1.4m)', 'LOFAR-LC (1.4m)', 'LOFAR-HC (1.4m)', 'LWA1 (3m)', 'LWA1x2x1', 'LWA1x4x1.5', 'LWA-OV (3m)', 'LWA-OVx2x1', 'LWA-OVx4x1.5']
    telescope_blmax = NP.asarray([5*14.0, 7*14.0, 21*14.0, 95*14.0, 100.0, 200.0, 1.4e3, 1.4e3, 1.4e3, 1.4e3, 300.0, 1e3, 1e3, 3.5e3, 3.5e3, 100.0, 100.0, 150.0, 200.0, 200.0, 300.0])
    telescope_wl = NP.asarray([2.0, 2.0, 2.0, 2.0, 0.5, 0.5, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 6.0, 2.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0])
    telescope_uvmax = telescope_blmax / telescope_wl
    telescope_n_antennas = NP.asarray([19, 37, 331, 6769, 1280, 1024, 112, 240, 496, 1008, 112, 0.75*1000, 0.75*256e3, 24, 24*2, 256, 512, 1024, 256, 512, 1024])
    # station_area = NP.asarray([NP.pi*(14.0/2)**2, NP.pi*(14.0/2)**2, 4.0**2, NP.pi*(35.0/2)**2, 3.0**2, 3.0**2, 1.38**2])
    telescope_antenna_area = NP.asarray([NP.pi*(14.0/2)**2, NP.pi*(14.0/2)**2, NP.pi*(14.0/2)**2, NP.pi*(14.0/2)**2, 20*(100.0/256), NP.pi*(6.0/2)**2, 4.0**2, 4.0**2, 4.0**2, 4.0**2, 4.0**2, NP.pi*(35.0/2)**2, 1.38**2, NP.pi*(87.0/2)**2, NP.pi*(30.8/2)**2, 3.0**2, 3.0**2, 3.0**2, 3.0**2, 3.0**2, 3.0**2])
    n_grid_telescopes_fov_gridding = telescope_blmax**2 / telescope_antenna_area
    ncomp_MOFF_telescopes_fov_gridding = 4*n_grid_telescopes_fov_gridding * NP.log2(4*n_grid_telescopes_fov_gridding)
    n_grid_telescopes_all_sky_gridding = 4 * telescope_blmax**2 / telescope_wl**2
    ncomp_MOFF_telescopes_all_sky_gridding = 4*n_grid_telescopes_all_sky_gridding * NP.log2(4*n_grid_telescopes_all_sky_gridding)
    ncomp_FX_telescopes = telescope_n_antennas * (telescope_n_antennas - 1) / 2
    print(telescopes[-6])
    print(ncomp_MOFF_telescopes_fov_gridding[-6])
    print(ncomp_FX_telescopes[-6])
    # ncomp_FX_telescopes = telescope_n_antennas ** 2
#     mrkrs = ['H','H','H','H','^','.','s','s','s','s','s','+','+','D','*','x','x','x','x','x','x']
#     msize = [4, 6, 8, 10, 4, 4, 4, 6, 8, 10, 12, 8, 12, 4, 6, 8, 10, 12, 14, 16, 18]
#     mew = [4, 4, 4, 4, 4, 8, 4, 4, 8, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4]
#     mfc = ['none', 'none', 'none', 'none', 'none', 'black', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none']

    # mrkrs = ['s', 'x', '*', '<', '>', 'v', '^', 'o', '.', '+', '+', '+', 'D', 'D', 'D']
    # msize = [4, 8, 8, 4, 4, 4, 4, 8, 4, 12, 14, 16, 4, 6, 8]
    # mew = [4, 4, 2, 4, 4, 4, 4, 4, 8, 4, 4, 4, 4, 4, 4]
    # mfc = ['none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'black', 'none', 'none', 'none', 'none', 'none', 'none']
    
#     fig = PLT.figure()
#     ax = fig.add_subplot(111)
#     ax.plot(10**NP.arange(11), 10**NP.arange(11), 'k:', lw=2)
#     for ti, telescope in enumerate(telescopes):
#         ax.plot(ncomp_MOFF_telescopes_fov_gridding[ti], ncomp_FX_telescopes[ti], mrkrs[ti], color='black', mfc=mfc[ti], ms=msize[ti], mew=mew[ti], label=telescope)
#     lgnd = ax.legend(loc='upper left', frameon=True, fontsize=10)
#     ax.set_xlim(0.1*ncomp_MOFF_telescopes_fov_gridding.min(), 10*ncomp_MOFF_telescopes_fov_gridding.max())
#     ax.set_ylim(0.1*ncomp_FX_telescopes.min(), 10*ncomp_FX_telescopes.max())
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.set_aspect('equal')
#     ax.set_xlabel('MOFF computations', fontsize=14, weight='medium')
#     ax.set_ylabel('FX computations', fontsize=14, weight='medium')
#     ax.xaxis.tick_top()
#     ax.yaxis.tick_right()
#     ax.xaxis.set_label_position('top')
#     ax.yaxis.set_label_position('right')    

#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_fov_gridding.png', bbox_inches=0)
#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_fov_gridding.eps', bbox_inches=0)    

#     fig = PLT.figure()
#     ax = fig.add_subplot(111)
#     ax.plot(10**NP.arange(11), 10**NP.arange(11), 'k:', lw=2)
#     for ti, telescope in enumerate(telescopes):
#         ax.plot(ncomp_MOFF_telescopes_all_sky_gridding[ti], ncomp_FX_telescopes[ti], mrkrs[ti], color='black', mfc=mfc[ti], ms=msize[ti], mew=mew[ti], label=telescope)
#     lgnd = ax.legend(loc='upper left', frameon=True, fontsize=10)
#     ax.set_xlim(0.1*ncomp_MOFF_telescopes_all_sky_gridding.min(), 10*ncomp_MOFF_telescopes_all_sky_gridding.max())
#     ax.set_ylim(0.1*ncomp_FX_telescopes.min(), 10*ncomp_FX_telescopes.max())
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.set_aspect('equal')
#     ax.set_xlabel('MOFF computations', fontsize=14, weight='medium')
#     ax.set_ylabel('FX computations', fontsize=14, weight='medium')
#     ax.xaxis.tick_top()
#     ax.yaxis.tick_right()
#     ax.xaxis.set_label_position('top')
#     ax.yaxis.set_label_position('right')    

#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_all-sky_gridding.png', bbox_inches=0)
#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_all-sky_gridding.eps', bbox_inches=0)    

    telescopes = ['HERA-19', 'HERA-37', 'HERA-331', 'HERA-6769', 'CHIME', 'HIRAX', 'MWA-112', 'MWA-240', 'MWA-496', 'MWA-1008', 'MWA-II-C', 'SKA1-LC', 'SKA1-LCD', 'LOFAR-LC', 'LOFAR-HC', 'LWA', 'LWA-x2x1', 'LWA-x4x2', 'LWA-OV', 'LWA-OVx2x1', 'LWA-OVx4x2']
    ind_HERA = NP.asarray([ti for ti,telescope in enumerate(telescopes) if telescope.split('-')[0] == 'HERA'])
    ind_MWA = NP.asarray([ti for ti,telescope in enumerate(telescopes) if telescope.split('-')[0] == 'MWA'])    
    ind_LWA = NP.asarray([ti for ti,telescope in enumerate(telescopes) if (telescope.split('-')[0] == 'LWA') or (telescope.split('-')[0] == 'LWA1')])
#     fig = PLT.figure()
#     ax = fig.add_subplot(111)
#     ax.plot(10**NP.arange(11), 10**NP.arange(11), 'k-', lw=2)
#     ax,plot(ncomp_MOFF_telescopes_fov_gridding[ind_HERA], ncomp_FX_telescopes[ind_HERA], 'k.', ls='--', lw=1)
#     ax,plot(ncomp_MOFF_telescopes_fov_gridding[ind_MWA], ncomp_FX_telescopes[ind_MWA], 'k.', ls=':', color='black', lw=2)    
#     # ax.plot(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[:2]], ncomp_FX_telescopes[ind_LWA[:2]], 'k--')
#     # ax.plot(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[2:]], ncomp_FX_telescopes[ind_LWA[2:]], 'k--')
#     ax.fill_betweenx(ncomp_FX_telescopes[ind_LWA[:3]], ncomp_MOFF_telescopes_fov_gridding[ind_LWA[:3]], ncomp_MOFF_telescopes_fov_gridding[ind_LWA[3:]], color='gray')
#     ax.annotate('LWA1', xy=(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[0]],ncomp_FX_telescopes[ind_LWA[0]]), xycoords='data', xytext=(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[0]],ncomp_FX_telescopes[ind_LWA[0]]/10), textcoords='data', arrowprops=dict(arrowstyle='->', color='black'), color='black', fontsize=9, horizontalalignment='center')
#     ax.annotate('LWA-OV', xy=(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[3]],ncomp_FX_telescopes[ind_LWA[3]]), xycoords='data', xytext=(ncomp_MOFF_telescopes_fov_gridding[ind_LWA[3]],ncomp_FX_telescopes[ind_LWA[3]]/10), textcoords='data', arrowprops=dict(arrowstyle='->', color='black'), color='black', fontsize=9, horizontalalignment='center')    
#     for ti, telescope in enumerate(telescopes):
#         if telescope.split('-')[0] != 'LWA':
#             ax.annotate(telescope, xy=(ncomp_MOFF_telescopes_fov_gridding[ti], ncomp_FX_telescopes[ti]), xycoords='data', horizontalalignment='center', verticalalignment='center', size=9)
#     ax.set_xlim(0.1*ncomp_MOFF_telescopes_fov_gridding.min(), 10*ncomp_MOFF_telescopes_fov_gridding.max())
#     ax.set_ylim(0.1*ncomp_FX_telescopes.min(), 10*ncomp_FX_telescopes.max())
#     ax.set_xscale('log')
#     ax.set_yscale('log')
#     # ax.set_aspect('equal')
#     ax.set_xlabel('MOFF computations', fontsize=14, weight='medium')
#     ax.set_ylabel('FX computations', fontsize=14, weight='medium')
#     ax.xaxis.tick_top()
#     ax.yaxis.tick_right()
#     ax.xaxis.set_label_position('top')
#     ax.yaxis.set_label_position('right')    

#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_fov_gridding_annotated.png', bbox_inches=0)
#     PLT.savefig('/data3/t_nithyanandan/project_MOFF/simulated/MWA/figures/MOFF_FX_computations_fov_gridding_annotated.eps', bbox_inches=0)    
