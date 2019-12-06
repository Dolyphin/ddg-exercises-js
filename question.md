# question 1
* 为什么
* Simplices 是如何选择上的？
answer: 在 function toggleSimplex(type, id)里面
## 那么 selected Simplex 被初始化成为一个什么？
answer: undefined 在 setup.js中
## 那么为什么不把它定义成 new Mesh()呢？geometry 也是定义成undefined


# Q2: isComplex函数的含义是什么？
什么样的集合能够被称为complex呢？看Slide和notes吧。cl(Complex)=Complex
## isPureComplex呢？
所有的complex都包含在度数为k的复形（complex)上

# Q3:关于buildVertexVector（）cannot convert "undefined" to int 的错误
solution1:看看 geometry中 matrix 是怎么用的。
*没有加algebra
看一看Solution中,DenseMatrix是如何应用的
SparseMatrix是用Triplet建立的。
