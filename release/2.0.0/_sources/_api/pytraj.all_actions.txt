.. _pytraj.commmon_actions:

pytraj.all_actions
==================

.. note:: Just need to import pytraj and use the method directly

.. code-block:: python
    
    import pytraj as pt
    pt.pca(traj)
    pt.surf(traj, mask='@CA')


.. container:: custom-index
    
    .. raw:: html
    
        <script type="text/javascript" src='../_static/cindex.js'></script>

- Most of the methods' names start with 'calc_' but use can ignore 'calc_'. For example::

   pt.calc_pca(...) is the same as pt.pca(...)


.. automodule:: pytraj.all_actions
    :members:
    :undoc-members:
    :show-inheritance:
    :inherited-members:
