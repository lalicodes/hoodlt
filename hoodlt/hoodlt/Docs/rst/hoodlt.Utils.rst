.. _ListUtilities:

Utilities and External Programs
===============================

The available utilities are classified into the following categories:

.. toctree::
   
   hoodlt.Utils.LatticeFunctions
   hoodlt.Utils.Geometry
   hoodlt.Utils.Visualization
   hoodlt.Utils.Misc

How to use those is best illustrated by Examples:

In this example, we compute the Voronoi cells for the :math:`\mbox{MgZn}_2` lattice, compute the gaussian Curvature and other geometric quantities and write the coordinates of the lattice in a format that can be used to be drawn by other programs. 

.. code-block:: python

        import numpy as np
        from hoodlt.OTM.MgZn2_OTM import OTMLatMgZn2Base12 
        # import the OTM MgZn2 class
        import hoodlt.Utils.voro_hoodlt as vh
        # import the voro_hoodlt utility`
        import hoodlt.Utils.Curvature  as cv
        # import the functions naming_convention, draw_facets
        from hoodlt.Utils.Lattice_Geometry import naming_convention, draw_facets

        # we are going to compute the Voronoi cell for base 0 (A NC) and 5 (B NC) 
        ind_fig = [(0, 'A'), (5, 'B')]
        save_result = False

        curv_metric = np.zeros(len(ind_fig))
        curv_topo = np.zeros(len(ind_fig))
        curv_r_topo = np.zeros_like(curv_topo)
        curv_r_metric = np.zeros_like(curv_metric)

        latc = OTMLatMgZn2Base12(3, 1.0, 0.9)
        # get gamma_crit
        gamma_crit = latc.gamma_crit[1]
        # generate an otm object at exactly gamma_crit
        # thus, we are computing the Voronoi cell exactly at gamma_c.
        latc = OTMLatMgZn2Base12(3, 1.0, gamma_crit)

        # this is the number of cells
        num_cells = latc.typ

        # get the voronoi object from the Voronoi class
        voro = vh.make_voro_pp(latc)
        # note that this will generate intermediate files in the directory

        for ind, ind_name in enumerate(ind_fig):

                ind_ws, ind_lab = ind_name
                
                # output the corresponding cell as a simple polyhedra   
                npa = voro.ws_simple(ind_ws)
                # output the corresponding cell as a PolyHedra object
                ws = voro.ws(ind_ws)
    
                n_vertices = len(ws.faces)
                n_faces = ws.vertices.shape[0]

                # here we print the properties of the Voronoi cell
                print('WS cell number', ind_ws)
                print('number of vertices', ws.num_vertices)
                print('number of edges', ws.num_edges)
                print('number of faces', ws.num_faces)
    
                # compute the Guassian curvature
                curv = cv.CurvPolyHedra(ws)
                curv_topo[ind] = curv.curv_topo
                curv_metric[ind] = curv.curv_metric
                curv_r_topo[ind] = curv.curv_reference_topo
                curv_r_metric[ind] = curv.curv_reference_metric
    
                # draw facets 
                if save_result:
                        draw_facets(voro, ind_ws)

                nm = naming_convention(voro, ind_ws)

                # save the points for subsequent use
                fid = open(nm + '.points', 'w')
                np.save(fid, voro.neighs[ind_ws])

                # save the radius for subsequent use
                fid = open(nm + '.radius', 'w')
                np.save(fid, voro.neighs_rad[ind_ws])

                # save the label, whether A or B for subsequent use
                fid = open(nm + '.label', 'w')
                np.save(fid, np.array([ind_lab]))     

        n_m = latc.typ
        print('\n')
        print('Print curvature')
        print(curv_topo)
        print(curv_metric)
        print('rel curvature topo', np.sum(n_m*curv_topo)/np.sum(n_m*curv_r_topo))
        print('rel curvature metric', np.sum(n_m*curv_metric)/np.sum(n_m*curv_r_metric))
