import streamlit as st
 
menu = st.sidebar.radio('***',
    (
    "Нелинейная краевая задача", 
    "Конечно-элементная аппроксимация",    
    "Итерационный метод", 
    "Одношаговый итерационный метод",    
    "Метод Ньютона",   
    "Упрощенный вариант",          
    )
)

if menu == "Нелинейная краевая задача":
    r"""
##### Нелинейная краевая задача

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
    """
if menu == "Конечно-элементная аппроксимация":
    r"""
##### Конечно-элементная аппроксимация

Пространства

* конечно-элементное пространство $V$ (кусрочно-гладкие фенкции)

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

**Вариационная постановка**

Найти  $y \in V_0$ такуе, что

$\begin{aligned}
 (k(y)  \operatorname{grad} y,  \operatorname{grad} v) - 
 (f(y,x), v) = 0 ,
 \quad \forall v \in V_0
\end{aligned}$

Общая запись

$\begin{aligned}
   F(y,v) = 0
\end{aligned}$

    """
  
if menu == "Итерационный метод":
    r"""
##### Итерационный метод

Итерационное уточнение приближенного решения

$\begin{aligned}
   y(x) = \sum_{i=1}^{n} y_i \varphi_i(x) 
\end{aligned}$

Две возможности

* итерационное нахождение коэффициентов $y_i , \ i = 1,2, \ldots, n$

* итерационное уточнение решения $y(x)$
    """

if menu == "Одношаговый итерационный метод":
    r"""
##### Одношаговый итерационный метод

**Задача**


$\begin{aligned}
   F(y,v) = 0
\end{aligned}$

для $y,v \in V_0$

**Итерационный метод**

$\begin{aligned}
   b_{k+1} \Big ( \frac{y^{k+1} - y^k}{\tau_{k+1}}, v \Big ) + F(y^k,v) = 0 ,
     \quad k =0,1,\ldots
\end{aligned}$

$b_{k+1}(y,v) - $ билинейная форма
    """

if menu == "Метод Ньютона":
    r"""
##### Метод Ньютона


**Модельная задача**

$\begin{aligned}
 F(y,v) = (k(y)  \operatorname{grad} y,  \operatorname{grad} v) -  (f(y,x), v)
\end{aligned}$

**Линеаризация по решению ($y$)**

$\begin{aligned}
 F(y^{k+1},v) & = (k(y^{k+1})  \operatorname{grad} y^{k+1},  \operatorname{grad} v) -  (f(y^{k+1},x), v) \\
 & \approx \big(k(y^k + k'(y^k)(y^{k+1}-y^k))  \operatorname{grad} (y^k + y^{k+1}-y^k),  \operatorname{grad} v \big) \\
 & -  \big (f(y^k + k'(y^k)(y^{k+1}-y^k),x), v) \\ 
 & \approx  \big(k(y^k)  \operatorname{grad} (y^{k+1}-y^k),  \operatorname{grad} v \big) +
 \big(k'(y^k)\operatorname{grad} y^k (y^{k+1}-y^k)),  \operatorname{grad} v \big)  + (k(y^k)  \operatorname{grad} y^k,  \operatorname{grad} v) \\
 & -
 \Big (\frac{\partial f}{\partial y}(y^k,x)(y^{k+1}-y^k),x), v \Big ) - \big (f(y^k,x), v)\\ 
\end{aligned}$

**Каноническая форма**

$\begin{aligned}
   b_{k+1} \Big ( \frac{y^{k+1} - y^k}{\tau_{k+1}}, v \Big ) + F(y^k,v) = 0 ,
     \quad k =0,1,\ldots
\end{aligned}$

В нашем случае

$\begin{aligned}
 & \tau_{k+1} = 1 \\
 & b_{k+1} (y,v) = (k(y^k)  \operatorname{grad} y,  \operatorname{grad} v) + (\operatorname{grad} k(y^k) y,  \operatorname{grad} v) - \Big (\frac{\partial f}{\partial y}(y^k,x) y, v \Big )
\end{aligned}$

    """

if menu == "Упрощенный вариант":
    r"""
##### Упрощенный вариант

**Проблема**

$ b_{k+1} (y,v) \ne  b_{k+1} (v,y) - $ несимметричная билинейная форма

* несимметричная матрица для коэффициентов $y_i , \ i = 1,2, \ldots, n$
* уравнение конвекции-диффузии

**Простейшая линеаризация**

$\begin{aligned}
b_{k+1} (y,v) = (k(y^k)  \operatorname{grad} y,  \operatorname{grad} v) - \Big (\frac{\partial f}{\partial y}(y^k,x) y, v \Big )
\end{aligned}$

    """
    

    
