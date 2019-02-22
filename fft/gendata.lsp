; Для генерации входных данных
(defparameter *N* (expt 2 20))
(defvar *set* 
  (make-array *N*
              :element-type '(complex single-float)
              :initial-element #c(0.0 0.0) ))

(defun garm(A omega i)
  (let ((arg (* omega i)))
    (* A (cis arg)) ))
(defun norm(ampl a s i)
  (complex (* ampl (exp (- (expt (/ (- i a) s) 2)))) 0))

(defun f (i)
;   (norm 1 (/ *N* 2) 10 i))
  (+ (garm 1 1 i) (garm 1 5 i) (garm 1 3 i) ))
; (cis i))


(loop :for i :from 0 :to (- *N* 1)
      :do (setf (elt *set* i) (f i) ))

(with-open-file (fd-out #p"./data/data.dat"
                        :direction :output
                        :if-exists :supersede)
  (format fd-out "# ~a~%" (length *set*))
  (loop :for cnum :across *set*
      :do (format fd-out "~f ~f~%" (realpart cnum) (imagpart cnum) )))
 
