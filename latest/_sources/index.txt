.. pytraj documentation master file, created by
   sphinx-quickstart on Mon Jun 22 19:09:20 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)

Welcome
=======

``pytraj`` is a Python front-end of the popular ``cpptraj`` package. Its aim is to expose
``cpptraj``'s funtions to Python's ecosystem. Enjoy.

.. raw:: html

   <div class="col-md-3">
   <h2>Overview</h2>

.. toctree::
   :maxdepth: 2

   overview

.. raw:: html

   </div>
   <div class="col-md-3">
   <h2>Documentation</h2>

.. toctree::
   :maxdepth: 1

   installation
   trajectory
   topology
   read_and_write
   tutorials/index
   analysis
   modify_traj
   atom_mask_selection
   trajectory_slice
   parallel
   whatsnew
   faq
   developer_guide
   misc
   conda
   cookbook
   trajectory_viewer
   api

.. raw:: html

   </div>
   <div class="col-md-3">
   <h2>Plot</h2>

.. image:: images/PCA_heart.png
   :target: tutorials/plot.html

|

**Indices and tables**

* `fork and contribute <https://github.com/Amber-MD/pytraj>`_
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. raw:: html

   </div>
   </div>
   </div>

.. raw:: html

   </div>
   <div class="col-md-3">
   <h3>Jupyter notebook</h3>

.. image:: http://jupyter.org/assets/jupyterpreview.png
   :target: http://jupyter.org/
   :height: 200


.. raw:: html

   </div>
   <div class="col-md-3">
   <h3><a href=trajectory_viewer.html> Trajectory visualization </a></h3>

    <script src="ngl.embedded.min.js">
    </script>
    
    <script>
    
      // adapted from NGL and MDAnalysis websites

      if( !Detector.webgl ) Detector.addGetWebGLMessage();
    
      NGL.mainScriptFilePath = "ngl.embedded.min.js";
    
      function onInit(){
          var stage = new NGL.Stage( "viewport" );
          stage.loadFile( "_static/1tsu.pdb" ).then( function ( o ) {
               o.addRepresentation( "cartoon" );
               o.addRepresentation( "licorice", { sele: "water" } );
               o.centerView();
          });

          stage.setTheme( "light" )
      
          window.addEventListener( "resize", function( event ){
             stage.handleResize();
          }, false );
      }
    
      document.addEventListener( "DOMContentLoaded", function() {
          NGL.init( onInit );
      } );
    
    </script>
    
    <div id="viewport" style="max-width:100%; height:400px;"></div>

.. raw:: html

   </div>
   <div class="col-md-3">
   <h3>Try pytraj online</h3>

.. image:: http://mybinder.org/images/logo.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. raw:: html

   </div>
