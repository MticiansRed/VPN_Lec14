import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Краевая задача", 
    "Конечно-элементное решение",    
    "Вариационная задача", 
    "Дискретная задача",    
    )
)
  
if menu == "Краевая задача":
    r"""
##### Краевая задача

**Уравнение**

$\begin{aligned} - 
 \operatorname{div} (k(u) \operatorname{grad} u ) = f(u,x),
 \quad x \in \Omega
\end{aligned}$

**Граничные условия**

$\begin{aligned}
  u(x) = 0,
  \quad x \in \partial\Omega
\end{aligned}$

**Нелинейность**

* коэффициент $k(u)$
* правая часть $f(u,x)$

    """    
if menu == "Конечно-элементное решение":
    r"""
##### Конечно-элементное решение

Пространства

* конечно-элементное пространство $V$ (кусочно-гладкие фенкции)

* подпространство $V_0 = \{ v \ | \ v \in V, \ v(x) = 0,
  \ x \in \partial\Omega \}$ 

Пробные функции в $V_0$ 

$\varphi_i(x) , \quad i = 1,2, \ldots, n$

$\begin{aligned}
 \varphi_i(x_j) = \left \{ \begin{array}{cc}
  1 ,  &  i = j \\
  0 ,  &  i \neq j \\
\end{array}
\right .
\end{aligned}$

Приближенное решение

$\begin{aligned}
   u(x) \approx y(x) = \sum_{i=1}^{n} y_i \varphi_i(x) ,
   \quad y_i = y(x_i),
   \quad i = 1,2, \ldots, n
\end{aligned}$

    """
    
if menu == "Вариационная задача":
    r"""
##### Вариационная задача

**Интегральное тождество**

$\begin{aligned}
 \int_{\Omega} k(y)  \operatorname{grad} y \operatorname{grad} v \, d x = 
 \int_{\Omega} f(y,x) v(x) \, d x
\end{aligned}$

для $y, v \in V_0$

**Вариационная постановка**

Найти  $y \in V_0$ такую, что

$\begin{aligned}
 (k(y)  \operatorname{grad} y,  \operatorname{grad} v) = 
 (f(y,x), v) ,
 \quad \forall v \in V_0
\end{aligned}$

Скалярное произведение

$\begin{aligned}
 (y,v) = \int_{\Omega} y(x) v(x) \, d x
\end{aligned}$
    """
    
if menu == "Дискретная задача":
    r"""
##### Дискретная задача

Конечно-элементное решение

$\begin{aligned}
   y(x) = \sum_{i=1}^{n} y_i \varphi_i(x) ,
   \quad i = 1,2, \ldots, n
\end{aligned}$

Задача для коэффициентов

$\begin{aligned}
 & \Big (k \big (\sum_{j=1}^{n} y_j \varphi_j \big ) \sum_{j=1}^{n} y_j \operatorname{grad} \varphi_j  \operatorname{grad} \varphi_i \Big ) = 
 \Big (f \big (\sum_{j=1}^{n} y_j \varphi_j, x \big ), \varphi_i \Big ) \\
& i = 1,2, \ldots, n
\end{aligned}$



Для вектора коэффициентов $y=\{y_i\}$

$\quad$ система нелинейных уравнений

$\begin{aligned}
 F(y) = 0 
\end{aligned}$

    """
    



























