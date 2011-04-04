(defpackage :base-test
  (:use :cl :base :lisp-unit))
(in-package :base-test)

(define-test cross
  (assert-equal (v 0 0 1) (cross (v 1) (v 0 1))))
