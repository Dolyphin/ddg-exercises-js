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

在direction-field-design中有

solution将indexElement()的参数去掉

采用 数组引用的方式vertexIndex[vertices]比较好
# Q4:如何实现Star
首先要知道 Star的含义 St(i)就是指i的邻域
* 分步骤进行实现
** 首先假设这只是一个点

## 实现了点的star，但是edge和face的star是什么呢？
查一下，notes
# Q5 Clousure是否有构造性的编法

# Q: isComplex如何定义？
A:每个simplex的face都在集合中

Q: halfedge的结构都指的是三角形网格吗？

# Q: 如何找到边向量？

how to find u and v

A:用geometry中的vector方法

cotan=u.dot(v)/u.cross(v).norm()

## Q: 0-form的定义是什么？

是表示成矩阵吗？

## Q:如何写科学的psudocode

## Q:如何写逻辑图

# Q: 如何将SICP SICM FDG 和DDG结合起来？

利用SICM和FDG做手写习题



## Q: form的几何意义？

# Q: form是行向量还是列向量，离散微分是左乘还是右乘



## Q:alphago的启示



## Q:学以致用的办法



## Q:为什么柯西，蒙日那样的人可以成为数学家。从工程师到数学家，数学也可以像音乐一样作为一种爱好。



#### Q： scheme如何查询函数体的定义，像python中的help函数一样？