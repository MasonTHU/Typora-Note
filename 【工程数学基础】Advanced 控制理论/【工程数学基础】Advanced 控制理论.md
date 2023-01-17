# 【工程数学基础】Advanced 控制理论



## 01 特征值eigenvalue与特征向量eigenvector

> DR_CAN
>
> [【工程数学基础】1_特征值与特征向量 bilibili](https://www.bilibili.com/video/BV1fx41137Zm/?spm_id_from=333.999.0.0&vd_source=c60eeae32275857559a9b94f7140b292)

### 1.1 定义

​		在线性代数中，对于一个给定的线性变换$\boldsymbol{A}$，它的特征向量$\boldsymbol{v}$经过这个线性变换的作用后，得到的新向量仍然与原来的$v$保持在同一条直线上，但其长度或方向也许会改变。即
$$
\boldsymbol{Av}=\boldsymbol{\lambda v}
$$
​		其中$\boldsymbol{\lambda}$为标量，即特征向量仍然与原来的$\boldsymbol{v}$保持在同一条直线上。但其长度或方向也许会改变。

---------------

### 1.2 背景

​		如果我们想求解如下状态方程组中$x_1, x_2$的值，怎么求？
$$
\begin{aligned}
&\frac{d x_1}{d t}=x_1+x_2 \\
&\frac{d x_2}{d t}=4 x_1-2 x_2
\end{aligned} \Rightarrow \frac{d}{d t}\left[\begin{array}{l}
x_1 \\
x_1
\end{array}\right]=\left[\begin{array}{cc}
1 & 1 \\
4 & -2
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right]
$$
​		反正我是不会求，因为$\frac{d x_1}{d t}和\frac{d x_2}{d t}$中都有$x_1和x_2$，它们耦合在一起很难去求解，所以我们需要利用特征值、特征向量以及对角矩阵的概念作为方法来计算。

> 耦合：一个系统里面的两个或以上的状态变量存在相互影响、相互关联的作用
>
> 2022.11.20 这里我完全是一知半解，差的还很远！
>
> 2022.11.21 做了后面几章后，理解了很多：之所以耦合的原因就是因为状态矩阵$\begin{bmatrix}1 & 1 \\4 & -2 \\ \end{bmatrix}$对$x_1和x_2$都有影响，如果能将其化作对角矩阵$\begin{bmatrix}\lambda_1 & 0 \\ 0 & \lambda_2 \\ \end{bmatrix}$就可以使$x_1和x_2$解耦 decouple

---

### 1.3 特征向量和特征值举例

​		比如$v_1$
$$
A=\begin{bmatrix}
1 & 1 \\4 & -2 \\ 
\end{bmatrix}
\qquad
v_1=\begin{bmatrix}
1 \\ 2 \\ 
\end{bmatrix} \\

Av_1=\begin{bmatrix}
1 & 1 \\4 & -2 \\ 
\end{bmatrix}
\begin{bmatrix} 
1 \\ 2 \\ 
\end{bmatrix}
=\begin{bmatrix} 
3 \\ 0 \\ 
\end{bmatrix}
$$
​		这里$Av_1和v_1$并不在同一条直线上，所以$v_1$并不是特征向量。

​		而对于如下$v_2$
$$
A=\begin{bmatrix}
1 & 1 \\4 & -2 \\ 
\end{bmatrix}
\qquad
v_2=\begin{bmatrix}
1 \\ 1 \\ 
\end{bmatrix} \\

Av_2=\begin{bmatrix}
1 & 1 \\4 & -2 \\ 
\end{bmatrix}
\begin{bmatrix} 
1 \\ 1 \\ 
\end{bmatrix}
=\begin{bmatrix} 
2 \\ 2 \\ 
\end{bmatrix} \\
\\ 
Av_2=2v_2
$$
​		所以$Av_2和v_2$在同一条直线上，其中$ {\color{red}v_2}$是特征向量，${\color{red}2}$是特征值。

### 1.4 求解特征向量及特征值

​		由定义
$$
\boldsymbol{Av}=\boldsymbol{\lambda v} \\
其中\lambda为特征值，v为特征向量 \\
$$
​		得
$$
Av-\lambda v=0 \\
(A-\lambda I)v=0 \\
其中I=\begin{bmatrix}
1 &  & 0 \\
 & \ddots &  \\
 0&  & 1 \\ 
\end{bmatrix}为单位对角矩阵，任何矩阵乘以单位对角矩阵仍为它本身 \\
则想要(A-\lambda I)v=0有非零解 \\
\Rightarrow |A-\lambda I|=0
$$
​		e.g.求解$A=\begin{bmatrix}1 & 1 \\4 & -2 \\ \end{bmatrix}$的特征向量及特征值
$$
A=\begin{bmatrix}
1 & 1 \\4 & -2 \\ 
\end{bmatrix}
\qquad
A-\lambda I=
\begin{bmatrix}
1-\lambda & 1 \\4 & -2-\lambda \\ 
\end{bmatrix} \\
|A-\lambda I|=0\Rightarrow
\begin{vmatrix}
1-\lambda & 1 \\4 & -2-\lambda \\ 
\end{vmatrix}=0 \\
(1-\lambda)(-2-\lambda)-1\times 4=0 \\
-2+2\lambda -\lambda +\lambda ^2 -4=0 \\
\lambda ^2+\lambda -6=0 \\
(\lambda-2)(\lambda+3)=0 \Rightarrow

\left\{
\begin{aligned}
&\lambda _1=2 \\
&\lambda _2=-3
\end{aligned}
\right. \\
------------------------------------
\\ 当\lambda_1 = 2时 \\
(A-\lambda I)v_1=
\begin{bmatrix}
-1 & 1 \\4 & -4 \\ 
\end{bmatrix} v_1 =0 \\
v_1对应\lambda_1的特征向量 \\
\begin{bmatrix}
-1 & 1 \\4 & -4 \\ 
\end{bmatrix} 
\begin{bmatrix}
v_{11} \\ v_{12} \\ 
\end{bmatrix}=0 \\

\left\{
\begin{aligned}
& -v_{11}+v_{12}=0 \\
& 4v_{11}-4v_{12}=0
\end{aligned}
\right. \Rightarrow v_{11}=v_{12} \\ 
则特征向量可以在v_{11}=v_{12}上任意取，比如\begin{bmatrix}
1 \\ 1 \\ 
\end{bmatrix} \\
------------------------------------ \\
当\lambda_2 = -3时 \\
(A-\lambda I)v_2=
\begin{bmatrix}
4 & 1 \\4 & 1 \\ 
\end{bmatrix} v_2 =0 \\
v_2对应\lambda_2的特征向量 \\
\begin{bmatrix}
4 & 1 \\4 & 1 \\ 
\end{bmatrix} 
\begin{bmatrix}
v_{21} \\ v_{22} \\ 
\end{bmatrix}=0 \\

\left\{
\begin{aligned}
& 4v_{21}+v_{22}=0 \\
& 4v_{21}+v_{22}=0
\end{aligned}
\right. \Rightarrow 4v_{21}+v_{22}=0 \Rightarrow 可取
\left\{
\begin{aligned}
& v_{21}=1 \\
& v_{22}=-4
\end{aligned}
\right. \\
特征向量为\begin{bmatrix}
1 \\ -4 \\ 
\end{bmatrix} \\
$$

### 1.5 应用

​		在背景中我们说明了想要解决的问题，求解微分方程组中$x_1和x_2$的值
$$
\begin{aligned}
&\frac{d x_1}{d t}=x_1+x_2 \\
&\frac{d x_2}{d t}=4 x_1-2 x_2
\end{aligned} \Rightarrow \frac{d}{d t}\left[\begin{array}{l}
x_1 \\
x_1
\end{array}\right]=\left[\begin{array}{cc}
1 & 1 \\
4 & -2
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right]
$$
​		下面我们看看怎么求。
$$
\frac{d}{d t}\left[\begin{array}{l}
x_1 \\
x_1
\end{array}\right]=\left[\begin{array}{cc}
1 & 1 \\
4 & -2
\end{array}\right]\left[\begin{array}{l}
x_1 \\
x_2
\end{array}\right] \\
\Rightarrow \frac{d x}{d t} = Ax \\
我们这里首先引入过渡矩阵P \\
P=[v_1, v_2] \qquad P:coordinate \ transformation \ matrix \ 过渡矩阵 \\
非常重要的一点：P是A的特征向量 \\
------------------------------------ \\
AP=A[v_1, v_2]= A
\begin{bmatrix}
v_{11} & v_{21} \\ v_{12} & v_{22} \\ 
\end{bmatrix} =

\begin{bmatrix}
A	\begin{bmatrix}
	v_{11} \\
	v_{12}
	\end{bmatrix}
A	\begin{bmatrix}
	v_{21} \\
	v_{22}
	\end{bmatrix} 
\end{bmatrix} \\

\because Av_1 = \lambda_1 v_1 \qquad Av_2 = \lambda_2 v_2  \\
\therefore
\begin{bmatrix}
A	\begin{bmatrix}
	v_{11} \\
	v_{12}
	\end{bmatrix}
A	\begin{bmatrix}
	v_{21} \\
	v_{22}
	\end{bmatrix} 
\end{bmatrix} \\ =

\begin{bmatrix}
\lambda_1	\begin{bmatrix}
	v_{11} \\
	v_{12}
	\end{bmatrix}
\lambda_2	\begin{bmatrix}
	v_{21} \\
	v_{22}
	\end{bmatrix}
\end{bmatrix} =

\begin{bmatrix}
\lambda_1 v_{11} & \lambda_2 v_{21} \\ \lambda_1 v_{12} & \lambda_2 v_{22} \\ 
\end{bmatrix} \\ =

\begin{bmatrix}
\ v_{11} & \ v_{21} \\ v_{12} & v_{22} \\ 
\end{bmatrix}

\begin{bmatrix}
\ \lambda_1 & \ 0 \\ 0 & \lambda_2 \\ 
\end{bmatrix}=P\Lambda(\Lambda是对角矩阵) \\

\Rightarrow AP=P\Lambda \\
------------------------------------ \\
同时对等式两边\times P^{-1} \\
\Rightarrow P^{-1}AP=P^{-1}P\Lambda=\Lambda \\
\Rightarrow P^{-1}AP=\Lambda \\
这样，我们通过过渡矩阵P将原矩阵A对角化
$$
​		回到$\frac{dx}{d t}=Ax$上
$$
\frac{dx}{d t}=Ax \\
\dot{x}=Ax \\
引入新的状态变量y \qquad 令x=Py \qquad 则\dot{x}=P\dot{y} \\
\therefore
\dot{x}=Ax \Rightarrow P\dot{y}=APy \\
两边同时\times P^{-1} \\
P^{-1}P\dot{y}=\underline{P^{-1}AP}y \\
\Rightarrow \dot{y}=\Lambda y=

\begin{bmatrix}
2 & 10 \\ 0 & -3 \\ 
\end{bmatrix}y \\

\Rightarrow 
\left\{
\begin{aligned}
& \dot{y_1}=2y_1 \\
& \dot{y_2}=-3y_2
\end{aligned}
\right.
\Rightarrow
\left\{
\begin{aligned}
& y_1=c_1e^{2t} \\
& y_2=c_2e^{-3t}
\end{aligned}
\right.\quad c_1,c_2为常数 \\

x=py=
\begin{bmatrix}
1 & 1 \\ 1 & -4 \\ 
\end{bmatrix}
\begin{bmatrix}
c_1e^{2t} \\ c_2e^{-3t} \\ 
\end{bmatrix}
=
\begin{bmatrix}
c_1e^{2t}+c_2e^{-3t} \\ c_1e^{2t}-4c_2e^{-3t} \\ 
\end{bmatrix}
$$

### 1.6 原文

<img src="【工程数学基础】Advanced 控制理论/image-20221120202750195.png" alt="image-20221120202750195" style="zoom: 50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120210937274.png" alt="image-20221120210937274" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120211126773.png" alt="image-20221120211126773" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120211412466.png" alt="image-20221120211412466" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120211820952.png" alt="image-20221120211820952" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120212005196.png" alt="image-20221120212005196" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221120212103862.png" alt="image-20221120212103862" style="zoom:50%;" />

## 02 卷积的拉普拉斯变换 Laplace Transform of Convolution 数学证明

> DR_CAN
>
> [【工程数学基础】4 卷积的拉普拉斯变换 Laplace Transform of Convolution 数学证明_ bilibili](https://www.bilibili.com/video/BV1fs411p7zD/?spm_id_from=333.788.recommend_more_video.-1&vd_source=c60eeae32275857559a9b94f7140b292)

### 2.1 定义

​		2.1.1 拉普拉斯变换 Laplace transform

​		
$$
X_{(s)}=\mathcal{L}\left[x_{(t)}\right]=\int_0^{\infty} x_{(t)} e^{-s t} dt
$$

> ​	拉普拉斯变换的积分下限为$-\infty$，但是在工程控制中，不需要研究时间0点以前的事情。

​		2.1.2 卷积 Convolution
$$
x(t) * h(t)=\int_0^t x(\tau) h(t-\tau) d\tau
$$

> ​    $*$  代表卷积运算符。卷积为线性时不变系统中，输入输出之间的关系，即系统的输入会对未来一段时间之内的系统输出产生影响

### 2.2 背景

​		对于一个线性时不变系统，输入是$x_(t)$，经过一个系统system后，输出为$y_(t)$。

<img src="【工程数学基础】Advanced 控制理论/image-20221121184604477.png" alt="image-20221121184604477" style="zoom:50%;" />

​		对于动态系统，我们可以用卷积描述输入输出关系，即$y_{(t)}=x_{(t)}*h_{(t)}=\int_0^t x(\tau) h(t-\tau) d \tau $。卷积在03章详细讲解。

​		而对于这样一个微分方程，显然不容易求解，就需要用到拉普拉斯变换了（把时域分析转化为复数域上分析）。

### 2.3 卷积的拉普拉斯变换证明

​		证明  $\mathcal{L}\left[x{(t)} * h{(t)}\right]=X{(s)} H{(s)}$
$$
\mathcal{L}\left[x_{(t)} * h_{(t)}\right] = 
\int_0^{\infty} 
	x_{(t)}*h_{(t)}e^{-st}
dt
\\ =
\int_0^{\infty} 
	\int_0^t 
		x(\tau) h(t-\tau) 
	d\tau e^{-st}
dt
\\
$$

> ​	$这里的d\tau和dt，就相当于二重积分\iint_D F{(t,\tau)}dA的积分区域dA=d\tau dt$
>
> ​	$d\tau积分$从0~t

<img src="【工程数学基础】Advanced 控制理论/image-20221121192513604.png" alt="image-20221121192513604" style="zoom:50%;" />

> ​	$dt积分$从0~$\infty$

<img src="【工程数学基础】Advanced 控制理论/image-20221121192721987.png" alt="image-20221121192721987" style="zoom:50%;" />

> ​	改变积分顺序

<img src="【工程数学基础】Advanced 控制理论/image-20221121200910447.png" alt="image-20221121200910447" style="zoom:50%;" />
$$
\int_0^{\infty} 
	\int_0^t 
		x(\tau) h(t-\tau) 
	d\tau e^{-st}
dt
\\=
\int_0^{\infty} 
	\int_\tau^\infty 
		\underline {x(\tau) h(t-\tau) 
	 e^{-st}
dt}d\tau \\
\\
let \quad t-\tau=u \\ 
then \quad t=u+\tau \\
dt = du + \cancel{d\tau=0} \\
t \in[\tau, \infty) \Rightarrow u=t-\tau \in[0, \infty) \\
x(\tau) h(t-\tau)e^{-st}dt \\
\\=
\int_0^{\infty} 
	\int_0^\infty 
	 x(\tau) h(u)e^{-s(u+\tau)}
	dud
\tau \\
=
\underline{
\int_0^{\infty}
	x(\tau)e^{-s\tau}
d\tau}
\quad \underline{
\int_0^{\infty}
	h(u)e^{-su}
du} \\
=X(s)H(s) \\
\\
\Rightarrow \mathcal{L}[X(t) * h(t)]=\mathcal{L}\left[X{(t)}\right] \mathcal{L}[H(t)]=X(s)H(s) \qquad Q.E.D.
$$


### 2.4 原文

<img src="【工程数学基础】Advanced 控制理论/image-20221121181131608.png" alt="image-20221121181131608" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221121182246962.png" alt="image-20221121182246962" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221121181707263.png" alt="image-20221121181707263" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221121181824151.png" alt="image-20221121181824151" style="zoom:50%;" />

##  03 线性时不变系统的冲激响应与卷积 

> DR_CAN
>
> [【工程数学基础】3 变声的基础原理 理解卷积的含义 线性时不变系统的冲激响应与卷积 哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1cs411W74f/?spm_id_from=333.999.0.0&vd_source=c60eeae32275857559a9b94f7140b292)

### 3.1 定义

- 线性时不变系统：LTI（linear time invariant) System，只是时间维度上的平移，时间本身并不是参数，不对结果产生影响，不以时间为转移。

​		定义运算符(operator)：$$o\left\{\cdot\right\}$$，使**input** $$o\left\{f(t)\right\}=x(t)$$ **output**

> ​	3.1.1 线性 Linear

$$
o\left\{f_1(t)+f_2(t)\right\}=x_1(t)+x_2(t) \\
o\left\{af(t)\right\}=ax(t) \\
\Rightarrow o\left\{a_1 f_1(t)+a_2 f_2(t)\right\}=a_1 x_1(t)+a_2 x_2(t) \\
\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad\qquad 叠加原理(Superposition\ Principal) \\
-------------------------------------------
$$

> ​	3.1.2 时不变 Time Invariant

$$
o\left\{f(t)\right\}=x(t) \Rightarrow o\left\{f(t-\tau)\right\}=x(t-\tau)
$$



- 冲激响应：Impulse Response
- 卷积：Convolution
- $x_{(t)}*h_{(t)}=\int_0^t x(\tau) h(t-\tau) d \tau $

### 3.2 背景

> [从“卷积”、到“图像卷积操作”、再到“卷积神经网络”，“卷积”意义的3次改变_哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1VV411478E/?spm_id_from=333.880.my_history.page.click&vd_source=c60eeae32275857559a9b94f7140b292)
>
> <img src="【工程数学基础】Advanced 控制理论/image-20221202224044062.png" alt="image-20221202224044062" style="zoom:50%;" />	
>
> ​		举个栗子：一个人在 **x** 时刻吃的东西，请问 **t** 时刻胃里剩下多少食物？
>
> ​		我们知道，吃东西在给胃里不断增加食物的同时，也会消化食物，这是一个动态的过程，求某一个时间点的变量是很困难的，但不是说没有办法，我们需要用两个函数来分别表达某一时刻的摄入量和消化量。以上图为例，f(t)代表t时刻胃里摄入的食物量，比如x=12点吃进去的食物就有f(12)=10克。g(t)代表消化率，代表t-x时刻后剩余多少食物，比如t=14点，也就是14-12=2h后，剩余g(14-12)=25%，则此时剩下的食物量就是$$f(12)*g(14-12)=10*0.25=2.5g$$。

> ​		对于一个线性时不变系统，每一个时间点都可能会发生一件有着同样规律的事件(例子中的g(t)，不以时间为转移)，而这个事件会对未来的某一时间点产生影响。那请问，众多事件的集合对同一时间点产生的影响又该如何获得呢？这个问题是整篇的核心。
>
> ​		从图中很好理解，每一个 x 对 t 的影响都是$$f(x)g(t-x)$$，我们只需要把0~t时间的所有x积分就可以得到t时刻之前所有的事件对这一点的影响，也就是$$\int_{0}^{t}f(x)g(t-x)dx$$。这个公式，我们称为卷积。

### 3.3 卷积的证明

​		卷积图像上怎么理解，卷积到底卷在了哪里，这个在b站上有很多的直观概念，我还是更喜欢这个，看着就“卷”。

<img src="【工程数学基础】Advanced 控制理论/image-20221203111521015.png" alt="image-20221203111521015" style="zoom:50%;" />

​		下面数学上证明为什么一个**输入函数**“卷”上一个**传递函数**就变成**输出函数**了：

​		首先是一个典型的二阶线性时不变系统，二阶弹簧阻尼振子模型：

<img src="【工程数学基础】Advanced 控制理论/image-20221203115620655.png" alt="image-20221203115620655" style="zoom:50%;" />

​		它的输入输出关系如下：

<img src="【工程数学基础】Advanced 控制理论/image-20221203115706727.png" alt="image-20221203115706727" style="zoom:50%;" />

​		在某一时刻给一个冲击力或者说激励(图中红色矩形)，之后弹簧就会按照(红色曲线)振动或者说响应。如果看一段时间内多次激励带来的影响，根据叠加原理，就相当于下图每一点之和：

<img src="【工程数学基础】Advanced 控制理论/image-20221203120357749.png" alt="image-20221203120357749" style="zoom:50%;" />

​		如果红色蓝色黑色的冲击足够小，矩形的短边$$\Delta T$$足够小，那加和就会转变为积分，$$\sum \Rightarrow \int$$。

​		想要解析一段时间内的多次激励和最后响应之间的关系，我们需要刻画一个非常理想的激励源，让弹簧震动的规律更清晰，从微积分的角度我们希望每一份的宽度无限小，但是又有很大的能量。

​		这里我们用Unit Impulse 单位冲击函数(又名Dirac Delta 狄拉克函数)来刻画这个冲击：
$$
Dirac Delta:
\begin{cases}
\delta(t)=0 & t\neq0 \ \  宽度为0 \\ \\
\int_{-\infty}^{\infty}\delta(t)dt=1 & \ \ \  \qquad面积为1 \\
\end{cases}
$$
​		如果理解不了可以用下面辅助函数：
$$
\delta_\Delta(t)=
\begin{cases}
\frac{1}{\Delta(T)} & 0<t<\Delta T \\
\quad 0 & else \\
\end{cases}
$$
<img src="【工程数学基础】Advanced 控制理论/image-20221203123315115.png" alt="image-20221203123315115" style="zoom:50%;" />
$$
\lim_{\Delta T \rightarrow 0 }\delta_\Delta(t)=\delta(t)
$$
​		对于这个理想的冲击，我们可以看到其响应有以下对应关系，既然是单位冲击响应，那我们很容易将其乘以系数A来获得不同比例的激励源。

<img src="【工程数学基础】Advanced 控制理论/image-20221203125218304.png" alt="image-20221203125218304" style="zoom:50%;" />

​		有了A，我们就可以将其和前面提到的更复杂的激励源联合起来

<img src="【工程数学基础】Advanced 控制理论/image-20221203125629124.png" alt="image-20221203125629124" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221203125751523.png" alt="image-20221203125751523" style="zoom:50%;" />

​		我们再来想一下为什么需要一个单位冲击，因为单位冲击的响应更容易刻画，用这个冲击去积分也就更容易表征各种其他的冲击。

​		根据叠加原理：
$$
When \quad t=i\Delta T \\
x(t)=\sum_{i=0}^{i}\Delta T f(i\Delta T)h_\Delta(t-i\Delta T) \\
When \quad \Delta T \rightarrow 0 \\
Because \quad \lim_{\Delta T \rightarrow 0 }\delta_\Delta(t)=\delta(t) \Rightarrow 
		\lim_{\Delta T \rightarrow 0 }h_\Delta(t)=h(t) \\
Then \quad \lim_{\Delta T \rightarrow 0 }\sum_{i=0}^{i}\Delta T f(i\Delta T)h_\Delta(t-i\Delta T)=
\lim_{\Delta T \rightarrow 0 }\sum_{i=0}^{i}\Delta T f(i\Delta T)h(t-i\Delta T) \\
Let\quad \Delta T = d\tau \quad i \Delta T=\tau \\
\Rightarrow x(t)=\int_{0}^{t}f(\tau)h(t-\tau)d\tau=f(t)*h(t) \qquad Q.E.D.
$$

### 3.4 原文

<img src="【工程数学基础】Advanced 控制理论/image-20221202194309404.png" alt="image-20221202194309404" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221202194518875.png" alt="image-20221202194518875" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221202200101899.png" alt="image-20221202200101899" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221202200457938.png" alt="image-20221202200457938" style="zoom:50%;" />

<img src="【工程数学基础】Advanced 控制理论/image-20221202200754994.png" alt="image-20221202200754994" style="zoom:50%;" />

## 04 欧拉公式的证明

> DR_CAN
>
> [【工程数学基础】5_如何证明宇宙第一美公式？？—欧拉公式证明_ 哔哩哔哩_bilibili](https://www.bilibili.com/video/BV1Ts411u7jA/?spm_id_from=333.999.0.0&vd_source=c60eeae32275857559a9b94f7140b292)

### 4.1 定义

<img src="【工程数学基础】Advanced 控制理论/image-20221203174305141.png" alt="image-20221203174305141" style="zoom:50%;" />

​		欧拉(Euler)公式：
$$
e^{i\theta}=cos\theta+isin\theta
$$
​		它还有一个更大名鼎鼎更有魅力的版本：
$$
e^{i\pi}+1=0
$$
​		一个简短的公式，包含自然底数e，虚数i，圆周率$$\pi$$。

​		By the way，e是Euler的首字母，但却不是欧拉发现的，是伯努利家族的人通过金融复利计算发现的。

### 4.2 背景

​		背景在此时此刻我无法写明，只能说在复变在信号处理及很多工控领域都有运用。

​		我们仔细观察$$cos\theta+isin\theta$$，和坐标轴是能扯上关系的，这就够了。

### 4.3 证明

​		We have
$$
i^2=-1
$$
​		Let
$$
f(\theta)=\frac{e^{i\theta}}{cos\theta+isin\theta}
$$
​		我们对上式求导
$$
\begin{aligned}
f^{\prime}(\theta) & =\frac{i e^{i \theta}(\cos \theta+i \sin \theta)-e^{i \theta}(-\sin \theta+i \cos \theta)}{(\cos \theta+i \sin \theta)^2} \\
& =\frac{i e^{i \theta} \cos \theta-e^{i \theta} \sin \theta+e^{i \theta} \sin \theta-e^{i \theta} i \cos \theta}{(\cos \theta+i \sin \theta)^2} \\
& =0
\end{aligned}
$$
​		Therefore
$$
f(\theta)\ is\ a\ constant
$$
​		We assume that
$$
\theta=0 \\
f(\theta)=f(0)= \frac{e^{i0}}{cos0+isin0}=1
$$
​		Hence 
$$
\frac{e^{i\theta}}{cos\theta+isin\theta}=1 \\
\Rightarrow e^{i\theta}=cos\theta+isin\theta \qquad Q.E.D.
$$

## 05 复变函数及欧拉公式

> DR_CAN
>
> [【工程数学基础】6_SinX=2? 复变函数 欧拉公式_哔哩哔哩_ bilibili](https://www.bilibili.com/video/BV1TW411z77n/?spm_id_from=333.788.recommend_more_video.-1&vd_source=c60eeae32275857559a9b94f7140b292)

### 5.1 定义

​		Complex 复数
$$
z=a+bi \quad i=\sqrt{-1}
$$
<img src="【工程数学基础】Advanced 控制理论/image-20221203185237413.png" alt="image-20221203185237413" style="zoom:50%;" />

-----------------

​		Assume
$$
z_1=a_1+b_1i \\
z_2=a_2+b_2i
$$
​		If we want 
$$
z_1 = z_2
$$
​		Then need
$$
a_1=a_2 \\
b_1=b_2
$$
​		According to Euler's Formula
$$
e^{i\theta}=cos\theta+isin\theta
$$
​		let it times r
$$
re^{i\theta}=rcos\theta+risin\theta
$$
<img src="【工程数学基础】Advanced 控制理论/image-20221203191109231.png" alt="image-20221203191109231" style="zoom:50%;" />

​		So a complex can be represented by
$$
z=re^{i\theta} \\with \quad r=\sqrt{a^2+b^2} \quad \theta=arctan \frac{a}{b}
$$

### 5.2 背景

​		是否存在
$$
sinx=2
$$

### 5.3 证明

​		我们知道如果$$x\in R$$，则不可能，因为正弦函数$$sinx \in [-1,+1]$$。

​		所以我们需要在$$x \in C$$ 中找答案。
$$
z=a+bi \\
sinz=2=C \quad C>1
$$

------------

​		利用欧拉公式，我们可以得到
$$
e^{iz}=cosz+isinz \tag{1}
$$
​		We let
$$
z=-z
$$
​		Then
$$
e^{-iz}=cosz-isinz \tag{2}
$$
​		We let
$$
(1)-(2) \Rightarrow e^{iz}-e^{-iz}=2isinz \\
sinz = \frac{e^{iz}-e^{-iz}}{2i}=C
$$
​		
$$
\because z=a+bi \\
\therefore sinz = \frac{e^{i(a+bi)}-e^{-i(a+bi)}}{2i}=\frac{e^{ai}e^{-b}-e^{-ai}e^{b}}{2i} \tag{3}
$$
​		再用一次欧拉公式
$$
e^{ai}=cosa+isina \\
e^{-ai}=cosa-isina
$$
​		则(3)式
$$
=\frac{e^{-b}(\cos{a}+i\sin{a})-e^{b}(\cos{a}-i\sin{a})}{2i} = \frac{(e^{-b}-e^{b})\cos{a}+(e^{-b}-e^{b})i\sin{a}}{2i} \\
=\frac{1}{2}(e^{-b}-e^{b})\sin{a} + [\frac{1}{2}(e^{b}-e^{-b})\cos{a}] i \\
=A+Bi=C \\
=C+0i
$$
​		Therefore

​		
$$
\left\{
\begin{aligned}
\frac{1}{2}(e^{-b}-e^{b})\sin{a} & = C \\
\\
\frac{1}{2}(e^{b}-e^{-b})\cos{a} & = 0 \\
\end{aligned}
\right.
$$
​		If
$$
\frac{1}{2}(e^{b}-e^{-b})\cos{a} = 0
$$
​		Then
$$
\Rightarrow
\left\{
\begin{aligned}
e^b-e^{-b} &= 0 \\
\\
\cos a    &= 0
\end{aligned}
\right.

\Rightarrow
\left\{
\begin{aligned}
b &= 0 \\
\\
a &= \frac{\pi}{2}+2k\pi
\end{aligned}
\right.
$$
​		**(1)	If	$$b = 0$$**
$$
\frac{1}{2}(1+1)\sin a=C \\
\Rightarrow \sin a = C \\
$$
​				But C>1

​				Therefore

​		
$$
\sin a \neq C \Rightarrow b \neq0 \qquad No\ Solution
$$
​		**(2)	If 	$$a = \frac{\pi}{2}+2k\pi$$**
$$
\frac{1}{2}(e^{-b}-e^{b})=C \\
e^{-b}-e^{b}-2C=0
$$
​				两边同时times $$e^b$$
$$
e^b \times (e^{-b}-e^{b}-2C)=e^b \times 0 \\
1+e^{2b}-2Ce^b=0
$$
​				Let
$$
e^b=u
$$
​				Then
$$
1+u^2-2Cu=0 \\
\Rightarrow u = \frac{2C \pm \sqrt{4C^2-4}}{2}=C \pm \sqrt{C^2-1}
$$
​				Then

​				
$$
b=\ln(c \pm \sqrt{C^2-1})
$$

----

​		Finally, if
$$
\sin z = C \quad C>1
$$
​		Then
$$
z=\frac{\pi}{2}+2k\pi + \ln(C \pm \sqrt{C^2-1})i
$$
​		If


$$
C=2
$$
​		Then
$$
z=\frac{\pi}{2}+2k\pi + \ln(2\pm \sqrt{3})i \qquad Q.E.D.
$$
