(require :asdf)
(require :raytrace)
(load "/home/martin/.local/share/common-lisp/source/lisp-unit/lisp-unit.lisp")
(use-package :lisp-unit)

(load "base-test.lisp")
(run-tests)