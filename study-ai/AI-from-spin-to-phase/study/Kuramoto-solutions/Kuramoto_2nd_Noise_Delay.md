> ## 1차 및 2차 쿠라모토 모형

- 외부 잡음(Noise)

- 시간 지연(Time Delay)

---

### 1. 1차 쿠라모토 모형 (관성 없음)

1차 모형은 진동자(오실레이터)가 상호작용 시 지연이나 물리적 저항(관성) 없이 위상을 즉각적으로 조정한다고 가정하는 모델입니다.

**① 결정론적 1차 모형 (기본형: 외부 잡음 없음, 시간 지연 없음)**

외부 교란이 없는 순수 비선형 시스템으로, 개별 진동자의 위상 변화율은 자신의 고유 진동수와 다른 진동자들과의 위상차에 의해서만 결정됩니다.

$$ \frac{d\theta_i(t)}{dt} = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) $$

**② 확률적 1차 SDE 모형 (외부 잡음 추가)**

결정론적 모형의 위상 갱신 과정에 시장 잡음이나 외부 교란을 나타내는 확률적 노이즈(표준 위너 프로세스 등)를 가산 형태로 결합한 이토 확률미분방정식(Itô SDE)입니다.

$$ d\theta_i(t) = \left[ \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] dt + \sigma dW_i(t) $$

**(※ 문헌에 따라 가우시안 백색 잡음 $\xi_i(t)$를 사용하여 $\dot{\theta}_i = \dots + \xi_i(t)$ 형태로 표기하기도 합니다)**

---

### 2. 2차 쿠라모토 모형 (관성 및 감쇠 포함)

2차 모형은 진동자에 회전 관성을 의미하는 **질량($m$)** 과 저항 역할을 하는 **감쇠 계수($\alpha$)** 를 추가하여 뉴턴 역학적인 가속도 운동을 기술합니다.

**① 결정론적 2차 모형 (기본형: 외부 잡음 없음, 시간 지연 없음)**

히스테리시스(이력 현상)와 1차 불연속 상전이 등 다중 안정성을 유발하는 결정론적 상미분방정식입니다.

$$ m \ddot{\theta}_i(t) + \alpha \dot{\theta}_i(t) = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) $$

**② 확률적 2차 SDE 모형 (외부 잡음 추가)**

2차 지배 방정식의 가속력(Force) 단에 예측 불가능한 가우시안 백색 잡음을 결합한 2계 확률미분방정식(Langevin SDE)입니다.

$$ m \ddot{\theta}_i(t) + \alpha \dot{\theta}_i(t) = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) + \sigma \xi_i(t) $$

**③ 2차 시간 지연 + 확률적 모형 (시간 지연 및 외부 잡음 동시 포함)**

정보나 자극이 즉각적으로 전달되지 않고 물리적 시차 ($\tau$)를 두고 전달되는 '리드-래그(Lead-Lag) 동역학' 및 '전도 지연' 현상을 반영한 가장 복잡한 형태의 모델입니다. 진동자 $i$는 현재가 아닌 **과거 $\tau$ 시점의 다른 진동자 위상 ($\theta_{j}(t-\tau)$) ** 에 반응합니다.

$$ m_i \frac{d^2 \theta_i(t)}{dt^2} + \alpha_i \frac{d\theta_i(t)}{dt} = \omega_i + \frac{1}{N} \sum_{j=1}^{N} K_{ij} \sin\bigl(\theta_j(t - \tau) - \theta_i(t)\bigr) + \sigma_i \xi_i(t) $$

---

### 💡 주요 수식 기호 설명

- **$\theta_i(t)$**: 진동자 $i$의 현재 위상

- **$\omega_i$**: 진동자 $i$의 고유 진동수 (독립적인 본연의 회전/변동 속도)
- **$K$ 및 $A_{ij}$ (또는 $K_{ij}$)**: 결합 강도 및 네트워크 연결 가중치(인접 행렬)
- **$m_i$**: 관성을 나타내는 질량 항 (변화에 저항하여 기존 추세를 유지하려는 힘)
- **$\alpha_i$**: 감쇠 계수 (움직임을 억제하는 마찰/저항력)
- **$\sigma, \xi_i(t), dW_i(t)$**: 외부에서 가해지는 예측 불가능한 잡음(노이즈)의 강도 및 함수
- **$\tau$**: 상호작용의 전달에 소요되는 물리적 시차 (시간 지연)

---

> ### 외부 잡음이 있을 때 시간 지연 2차 Kuramoto 모델의 구조와 풀이 방법을 수식을 이용해서 쉽게 설명해줘

- 외부 잡음(Langevin Noise)과 시간 지연(Time Delay)이 존재하는 **2차 Kuramoto Stochastic Differential Equation (SDE) 모델** 의 수학적 구조와 이를 컴퓨터로 해결하는 수치적 풀이 방법을 수식을 곁들여 설명해 드리겠습니다.

---

## 1. 모델의 수학적 구조 (Equation Structure)

시간 지연이 있는 2차 Kuramoto SDE 모델은 각 진동자(오실레이터) $i$ ($i = 1, 2, \dots, N$)에 대해 다음과 같은 물리 방정식으로 정의됩니다.

$$m_i \frac{d^2 \theta_i(t)}{dt^2} + \alpha_i \frac{d\theta_i(t)}{dt} = \omega_i + \frac{1}{N} \sum_{j=1}^{N} K_{ij} \sin\bigl(\theta_j(t - \tau) - \theta_i(t)\bigr) + \sigma_i \xi_i(t)$$

각 항이 주가 예측 모델에서 갖는 물리적 의미는 다음과 같이 대입해 이해할 수 있습니다.

| 항 | 공식 | 물리적 의미 | 금융학적 해석 |
| :--- | :--- | :--- | :--- |
| **관성 항** | $m_i \frac{d^2 \theta_i}{dt^2}$ | 질량($m_i$)에 의한 움직임 가속 저항 | 주가가 직전의 상승/하락 추세를 유지하려는 **관성(Momentum)** |
| **점성 감쇠 항** | $\alpha_i \frac{d\theta_i}{dt}$ | 움직임의 속도를 억제하는 마찰력 | 주가가 과도하게 요동치지 않도록 브레이크를 밟아주는 **감쇠 및 안정 장치** |
| **자연 주파수** | $\omega_i$ | 고유의 회전 속도 | 시장 외력이 없을 때 주가가 본래 흘러가려는 **장기 고유 추세(Base Trend)** |
| **지연 결합 항** | $\frac{1}{N} \sum K_{ij} \sin(\theta_j(t-\tau) - \theta_i)$ | **$\tau$시간 전 과거** 진동자들의 위상과 상호작용 | 삼성전자 등 주도주가 변동한 후 **$\tau$일의 시차(Time Delay)**를 두고 주변 주로 충격이 전파되는 **Lead-Lag 동역학** |
| **외부 잡음 항** | $\sigma_i \xi_i(t)$ | 무작위로 타격을 주는 Langevin 노이즈 | 거시 뉴스, 돌발 호재/악재 등 예측 불가능한 **시장 불확실성(Random Walk 변동)** |

---

## 2. 수치적 풀이 방법 (Numerical Solution Method)

컴퓨터는 연속적인 미분 방정식을 직접 풀 수 없기 때문에, 아주 미세한 시간 보폭 $dt$ (예: `0.05`초) 단위로 쪼개어 미래 상태를 한 스텝씩 전진 계산하는 **수치적 적분(Numerical Integration)** 방식을 취합니다.

여기서는 **Runge-Kutta 4th-order (RK4)** 기법과 확률적 잡음(Stochastic Noise)을 결합해 해결합니다.

### 1단계: 2차 미분 방정식을 1차 연립 미분 방정식으로 변환

컴퓨터 처리를 쉽게 하기 위해 각속도(주가 변동 속도) $v_i(t) = \frac{d\theta_i(t)}{dt}$라는 새로운 변수를 도입하여, 2차 미분 방정식을 **두 개의 1차 연립 미분 방정식**으로 나눕니다.

1. $\frac{d\theta_i(t)}{dt} = v_i(t)$

2. $\frac{dv_i(t)}{dt} = \frac{1}{m_i} \left[ \omega_i - \alpha_i v_i(t) + \frac{1}{N} \sum_{j=1}^{N} K_{ij} \sin\bigl(\theta_j(t - \tau) - \theta_i(t)\bigr) \right] + \frac{\sigma_i}{m_i} \xi_i(t)$

---

### 2단계: 결정론적(Deterministic) 부분의 RK4 적분

잡음 $\xi_i(t)$를 제외한 순수 물리 도함수(기울기)를 각각 $f_\theta = v_i(t)$, $f_v = \text{Right Hand Side}$ 라고 하겠습니다.

RK4 기법은 현재 시점 $t$에서 다음 시점 $t+dt$로 가기 위해, 중간 시점의 기울기($k_1, k_2, k_3, k_4$)를 네 차례 밟아 평균을 내어 전진시킵니다.

1. **$k_1$ (현재 위치에서의 기울기)**:

   - $k_{1\theta} = v(t)$
   - $k_{1v} = f_v\bigl(\theta(t), v(t), \theta(t-\tau)\bigr)$
     
3. **$k_2$ (절반 보폭 $dt/2$만큼 전진한 위치의 기울기)**:

   - $k_{2\theta} = v(t) + 0.5 \cdot dt \cdot k_{1v}$
   - $k_{2v} = f_v\bigl(\theta(t) + 0.5 \cdot dt \cdot k_{1\theta}, v(t) + 0.5 \cdot dt \cdot k_{1v}, \theta(t-\tau)\bigr)$
     
4. **$k_3$ (다시 그 중간 위치에서 새로 조정한 기울기)**:

   - $k_{3\theta} = v(t) + 0.5 \cdot dt \cdot k_{2v}$
   - $k_{3v} = f_v\bigl(\theta(t) + 0.5 \cdot dt \cdot k_{2\theta}, v(t) + 0.5 \cdot dt \cdot k_{2v}, \theta(t-\tau)\bigr)$
     
5. **$k_4$ (전체 보폭 $dt$만큼 전진한 끝점에서의 기울기)**:

   - $k_{4\theta} = v(t) + dt \cdot k_{3v}$
   - $k_{4v} = f_v\bigl(\theta(t) + dt \cdot k_{3\theta}, v(t) + dt \cdot k_{3v}, \theta(t-\tau)\bigr)$

- **지연 항 처리**: 과거 위상 $\theta_j(t-\tau)$는 이 단계들 중에는 고정된 과거 기록이므로, 수치 연산 전에 메모리 이력 버퍼(History Buffer)에서 **선형 보간(Linear Interpolation)** 을 통해 값을 찾아내어 상수로 대입합니다.

---

### 3단계: 외부 Langevin 노이즈 항(Stochastic Term) 주입

위 단계들을 가중 평균하여 잡음이 없는 다음 단계 예측치 $(\theta_{new}, v_{new})$를 구합니다.

- $$\theta_{new} = \theta(t) + \frac{dt}{6} (k_{1\theta} + 2k_{2\theta} + 2k_{3\theta} + k_{4\theta})$$

- $$v_{new} = v(t) + \frac{dt}{6} (k_{1v} + 2k_{2v} + 2k_{3v} + k_{4v})$$

여기에 확률 변동성인 **Langevin 노이즈**를 주입합니다. 확률 적분 법칙상, 시간 간격 $dt$ 동안 브라운 운동의 무작위 증분은 $\sqrt{dt}$ 스케일을 따르므로 표준 정규분포 난수 $Z \sim \mathcal{N}(0, 1)$를 생성하여 최종적으로 다음 단계 값을 얻습니다.

- $$v(t + dt) = v_{new} + \frac{\sigma_i}{m_i} \sqrt{dt} \cdot Z_v$$
  
- $$\theta(t + dt) = \theta_{new} + 0.1 \cdot \sigma_i \sqrt{dt} \cdot Z_\theta$$

이 과정을 시뮬레이션 전체 일수만큼 루프를 돌며 계속 반복하는 것이 **시간 지연 2차 Kuramoto SDE 모델의 수치적 풀이 과정**입니다.

---
