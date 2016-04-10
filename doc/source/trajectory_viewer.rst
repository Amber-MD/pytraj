Trajectory Viewer
=================

.. note:: 

    Still work in progress. Please see source code in https://github.com/arose/nglview/

    or `try this online <http://mybinder.org/repo/hainm/notebook-pytraj/>`_

Requirement: jupyter notebook, nglview::

    # install notebook via conda.
    conda install jupyter notebook

    # install nglview
    git clone https://github.com/arose/nglview
    cd nglview
    python setup.py install

| 
.. raw:: html

    <script src="_static/ngl.embedded.min.js">
    </script>
    
    <script>
    
      if( !Detector.webgl ) Detector.addGetWebGLMessage();
    
      NGL.mainScriptFilePath = "_static/ngl.embedded.min.js";
    
      function onInit(){
          var stage = new NGL.Stage( "viewport" );
          stage.loadFile( "_static/1tsu.pdb", { defaultRepresentation: true } );
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
