# SlicerMitralValve

This extension contains a collection of tools for aiding in patient-specific mitral valve modelling.

## Biplane Registration

A scripted module that partially automates extracting 2 image planes from Philips bi-plane ultrasound and aligning them in 3D space. Allows for minimal user input to reach final registration.

## MVModeller (Deprecated)

Loadable module prototyping a manual segmentation/modelling process in Slicer. Replaced by MVSegmenter module implementing automated process.

## MVSegmenter

Scripted module implementing automatic mitral valve segmentation using ITK. Depends on HeartValveLib from the Slicer Heart extension.
