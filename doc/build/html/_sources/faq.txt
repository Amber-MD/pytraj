.. _faq:

Frequently Asked Questions
==========================

- Why do most of examples use **pytraj.iterload** instead of **pytraj.load**

    Because we encourage user to use out-of-core calculation for memory saving. Most of
    the time we don't need change trajecotry's coordinate, so it's OK to just use
    immutable trajectory.
