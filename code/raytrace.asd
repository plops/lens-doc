(asdf:defsystem raytrace
  :components
  ((:file "packages")
   (:file "base" :depends-on ("packages"))
   (:file "raster" :depends-on ("packages" "base"))
   (:file "raytrace" :depends-on ("packages" "base"))
   (:file "plot-macros" :depends-on ("packages" "base"))
   (:file "project" :depends-on ("packages" "base" "raytrace"))))
