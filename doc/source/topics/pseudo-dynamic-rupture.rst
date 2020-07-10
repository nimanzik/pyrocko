**************************************************************************
Pseudo Dynamic Rupture - *A stress-based self-similar finite source model*
**************************************************************************

Introduction
============

Physics-based, dynamic rupture models, which rely on few parameters only, are needed to allow a realistic forward modeling and to reduce the effort and non-uniqueness of the inversion. Our study develops a simplified, dynamic, semi-analytical rupture model suited for wavefield simulation and inversion. The rupture builds on the class of self-similar crack models. On one hand it is approximative as it neglects inertia and so far the details of frictional effects, and treats the rupture front growth in a simplified way.  On the other hand, it is complete as the slip fulfils the boundary condition on the broken plane for every instantaneous rupture front geometry and applied stress. 

.. contents :: Content
  :depth: 2

Source Implementation details
=============================

.. figure :: /static/pseudo-dynamic-flow.svg
    :align: center
    :width: 90%
    :alt: Flowchart calculation of the pseudo dynamic rupture.

    Flowchart illustrating the implementation of the ``PseudoDynamicRupture`` in Pyrocko-GF.

Pyrocko-GF source usage
=======================

The ``PseudoDynamicRupture`` can be used a any other source model in Pyrocko-GF for forward modelling synthetic waveforms or static surface displacements.

Traction and stress field parametrisation
-----------------------------------------

The rupture plane can be exposed to different stress/traction field models and abstract modifiers to these fields. These fields can be used independently or be combined into a composition model.

Available models:

    * ``UniformTractions``
    * ``HomogeneousTractions``

Abstract models:

    * ``RectangularTaper``
    * ``DepthTaper``

Forward calculation of waveforms and static displacement
--------------------------------------------------------

Details how the source model can be used in Pyrocko-GF can are layed out here: :doc:`pyrocko-gf`.

Plotting source parameters
==========================

Slip distribution
-----------------

Rupture evolution
-----------------

Moment release
--------------
