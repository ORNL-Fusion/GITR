# Release Notes

This first set of release notes details new features added to GITR.

## Automatic dependency handling

Dependencies are installed automatically if not found on host system.

## Reduced repo size

Pre-processing utilities and examples moved to different repositories.
Check GITR_legacy_<example_name> under ORNL/Fusion for examples.
Check GITR_processing under ORNL/Fusion for processing utilities

## New directory layout

New, flat directory structure implemented.

## Header files split into declaration/definition

Many modules have been split into header/implementation. Previously
these were all defined in headers.

## Disabled defunct unit tests
Unit tests that tested unused functionality have been removed.
