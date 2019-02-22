; тут пока меньше 200 строк, а разделять файлы -- сложно
; ещё утвержается, что передавать массивы между независимо компилируемыми
;+ файлами неэффективно
; а ещё получился почти фортран, что весьма показательно
; этот код наверное можно снабжать повестью
; 1) 2^10 -- 4s    (gcl, interp)
; 2) 2^10 -- 0.90s (sbcl, comp)
; 3) 2^10 -- 0.70s (sbcl, comp) -- убраны замыкания
; 4) 2^10 -- 0.20s (sbcl, comp) -- предпосчитаны W (лениво)
; 5) 2^12 -- 0.16s (sbcl, comp) -- убраны остатки функциональщины
; 6) 2^16 -- 72s (было 78) (sbcl, comp) -- мелкие докрутки
; 7) 2^16 -- 48s (sbcl, comp) -- полный предпосчёт W
; 8) пришло осознание, что исходный алгоритм работал за N^2 * log(N)
; 9) 2^20 -- 5s (sbcl, comp) 
; 10) больше не получится -- лисп пишет, что кончается куча..

; все оптимизировать!
(declaim (optimize (safety 0) (space 0) (speed 3) (debug 0)))

; куча констант, нужны именно тут для оптимизаций при компиляции
(deftype uint_32 () `(integer 0 ,(expt 2 32))  )
(deftype ctype   () '(complex single-float)    )
(deftype cartype () '(simple-array ctype (*))  )
(deftype fltype  () 'single-float              )
(deftype fltypepos () '(single-float 0.0)     )
(defconstant +c-zero+ (coerce #c(0 0) 'ctype)) (declaim (type ctype +c-zero+))
(defconstant +PI+     (coerce pi 'fltype) ) (declaim (type fltype +PI+))

(defvar *data*)      (declaim (type cartype *data*    ))
(defvar *res* )      (declaim (type cartype *res*     ))
(defvar *numpoints*) (declaim (type uint_32 *numpoints*))
(defvar *W* ) (declaim (type (simple-array ctype (* *)) *W*))
(defconstant +quadlim+ 2)(declaim (type fixnum +quadlim+)) 

; идея с длиной числа в битах утащена из maxima : fft-core.lisp
(declaim (ftype (function (uint_32) uint_32) up-pow-2) (inline up-pow-2))
(defun up-pow-2 (x)
  (let ((low (ash 1 (1- (integer-length x))) ))
    (if (= low x) low (ash low 1))))

(declaim (ftype (function (uint_32) uint_32) log-ind) (inline log-ind))
(defun log-ind(p) (- (integer-length p) 1))

(defmacro v (oper v1 v2)
  `(map 'cartype #',oper ,v1 ,v2))

; считаем все W
(defun calc-W (&key (direct T))
(let ((a (array-dimension *W* 0))
      (b (array-dimension *W* 1)))
  (declare (type uint_32 a) (type fixnum b))
  (dotimes (i a)
    (dotimes (j b)
      (declare (type uint_32 i) (type fixnum j))
      (let* ((q i)
             (p (ash 1 j))
             (arg (* (if direct -1 1) (/ (* 2 +PI+ q) p)) ))
        (declare (type uint_32 q p))
        (setf (aref *W* i j) (cis arg)) )))))

(declaim (ftype (function (cartype uint_32) cartype) quad-transform))
(defun quad-transform (x size)
  (let ((y (init-cmplx-arr size)))
    (declare (type cartype y))
    (loop :for k :of-type uint_32 :from 0 :to (1- size)
          :do (loop :for n :of-type uint_32 :from 0 :to (1- size)
                    :sum (* (aref x n) (aref *W* (* n k)(log-ind size)))
                    :into res :of-type ctype
                    :finally (setf (aref y k) (/ res (sqrt (the uint_32 size)))) ))
    y))


(declaim (ftype (function (cartype uint_32) cartype) fast-transform))
(defun fast-transform (x size)
  (if (> size +quadlim+)
    (let* ((x1     (subseq x 0              (floor size 2)) )
           (x2     (subseq x (floor size 2) size)           )
           (x1+x2  (v + x1 x2)                              )
           (x1-x2  (v - x1 x2)                              )
           (size/2 (floor size 2)                         )
           (y      (init-cmplx-arr size)                  ))
      (declare (type uint_32 size/2) (type cartype x1 x2 x1+x2 x1-x2 y))
      ; чётные
      (setf x1 (fast-transform x1+x2 size/2))
      (loop :for i :of-type uint_32 :from 0 :to (1- size/2)
            :do (setf (aref y (* 2 i)) (aref x1 i)) )
      ;нечётные
      (loop :for n :of-type uint_32 :from 0 :to (1- size/2)
            :do (setf (aref x1-x2 n)
                      (* (aref *W* n (log-ind size)) (aref x1-x2 n) )))
      (setf x2 (fast-transform x1-x2 size/2))
      (loop :for i :of-type uint_32 :from 0 :to (1- size/2)
            :do (setf (aref y (1+ (* 2 i))) (aref x2 i)) )
      y)
    (quad-transform x size)))

(declaim (ftype (function (cartype &key (:direct T)) cartype) fft))
(defun fft (x &key (direct T))
  (calc-W :direct direct)
  (fast-transform x (array-total-size x)))


(declaim (ftype (function (uint_32) cartype) init-cmplx-arr))
(defun init-cmplx-arr (size)
  (let ((arr (make-array size 
                         :element-type 'ctype
                         :initial-element +c-zero+)))
    arr))

(declaim (ftype (function (uint_32 uint_32) (simple-array ctype (* *)))
                          init-cmplx-mtx))
(defun init-cmplx-mtx (dim1 dim2)
  (let ((mtx (make-array (list dim1 dim2)
                         :element-type 'ctype
                         :initial-element +c-zero+)))
    mtx))



; (defun main()
(format t "Кажется, оно завелось~&")
 
(with-open-file (fd-in #p"./data/data.dat" :direction :input)
  (read-char fd-in); решётку
  (let* ((numpoints (read fd-in))
         (size      (up-pow-2 numpoints)))
    (declare (type uint_32 numpoints)
             (type uint_32 size))
    (setf *numpoints* numpoints)
    (setf *data* (init-cmplx-arr size))
    (setf *res* (init-cmplx-arr size))
    (setf *W* (init-cmplx-mtx 
                (the uint_32 (* +quadlim+ size)) 
                (1+ (log-ind size)) ))
    ; заполнение
    (loop :for i :of-type uint_32 :from 0 :to (1- *numpoints*)
          :do   ;поможем лиспу понять, как считывать
          (let ((buf (concatenate 'string "#c(" (read-line fd-in) ")" )))
            (setf (aref *data* i)  ;(read (make-string-input-stream buf)) 
              (coerce (read (make-string-input-stream buf))
                      (array-element-type *data*)) )))))

(format t "И даже прочитало данные~&")
(print (time (progn (setf *res* (fft *data*)) nil)))
(format t "~&Кажется, оно что-то посчитало~&")

(with-open-file (fd-out #p"./data/res.dat"
                        :direction :output
                        :if-exists :supersede
                        :if-does-not-exist :create)
  ; он вроде должен сам упасть если кончилось место
  (declare (type stream fd-out))
  (format fd-out "# ~a~%" *numpoints*)
  (loop :for i :of-type uint_32 :from 0 :to (1- *numpoints*)
        :do (let ((cnum (aref *res* i) ))
              (declare (type ctype cnum))
              (format fd-out "~f ~f~%" (realpart cnum) (imagpart cnum) ))))
(format t "~&И даже записало ответ~&")
;)

