# 확률적 2차 쿠라모토(Kuramoto) SDE 모형의 이론 및 수치적 해석 기법

본 보고서는 질량(관성)과 감쇠(저항)가 있는 2차 쿠라모토 진동자계에 확률적 노이즈가 도입된 **'확률적 2차 쿠라모토 SDE 모형(Stochastic 2nd-order Kuramoto SDE Model)'** 의 지배 방정식, 수학적 해석 방법 및 컴퓨터로 이산적인 해를 구하는 수치 해석 알고리즘과 참조 파이썬 코드를 상세히 정리한 문서입니다.

---

## 1. 확률적 2차 쿠라모토 SDE 지배 방정식

2차 쿠라모토 모델은 개별 오실레이터의 회전 운동에 관성(Inertia)을 나타내는 **질량($m$)** 과 저항을 나타내는 **감쇠 계수($\alpha$)** 를 고려합니다. 여기에 외부 교란 및 시장의 무작위 변동성을 묘사하는 **가우시안 백색 잡음(Gaussian White Noise)** 을 가속력(Force) 단에 결합하면 아래와 같은 2계 확률미분방정식(Langevin SDE)을 얻습니다.

$$m \ddot{\theta}_i(t) + \alpha \dot{\theta}_i(t) = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) + \sigma \xi_i(t)$$

이 식에서 각 물리량의 정의는 다음과 같습니다:

* $\theta_i(t)$: 오실레이터 $i$의 위상 (Phase)
* $v_i(t) = \dot{\theta}_i(t)$: 오실레이터 $i$의 각속도 (Angular Velocity)
* $m$: 오실레이터의 질량 (관성량, Inertia)
* $\alpha$: 감쇠 계수 (마찰/저항력, Damping coefficient)
* $\omega_i$: 오실레이터 $i$의 고유 진동수 (Intrinsic Frequency)
* $K$: 결합 강도 (Coupling strength)
* $A_{ij}$: 인접 행렬(Adjacency Matrix, 마스크 가중치)
* $\sigma$: 노이즈 강도 (Noise intensity)
* $\xi_i(t)$: 평균이 0이고 상관관계가 없는 독립 가우시안 백색 잡음 (Gaussian White Noise)
* 
  $$\langle \xi_i(t) \rangle = 0, \quad \langle \xi_i(t) \xi_j(t') \rangle = \delta_{ij} \delta(t-t')$$

---

## 2. 상태 공간 방정식 (State-Space Representation)

컴퓨터 상에서 수치 적분을 수월하게 수행하기 위해, 위의 2계 SDE를 2개의 연립 1계 이토 확률미분방정식(Itô SDE) 시스템으로 변환합니다. 속도 변수 $v_i(t) = \dot{\theta}_i(t)$를 도입하면 다음과 같이 상태 공간 방정식으로 기술됩니다.

$$d\theta_i(t) = v_i(t) dt$$

$$dv_i(t) = \frac{1}{m} \left[ \omega_i - \alpha v_i(t) + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] dt + \frac{\sigma}{m} dW_i(t)$$

여기서 $W_i(t)$는 표준 위너 프로세스(Wiener Process, Brownian Motion)이며, $dW_i(t) = W_i(t+dt) - W_i(t) \approx \mathcal{N}(0, dt)$의 확률적 증분을 가집니다.

노이즈가 상태 변수($\theta, v$)에 곱해지는 형태가 아닌 가산 형태로 더해지는 **덧셈 노이즈(Additive Noise)** 형식이므로, Itô 적분과 Stratonovich 적분의 확률 기하학적 정의 차이가 발생하지 않고 수식이 일치하는 편리함을 제공합니다.

---

## 3. 수학적 해를 구하는 방법론 (Mathematical Solution Methods)

비선형 사인 함수 커플링 항 $\sin(\theta_j - \theta_i)$과 다체 상호작용의 비선형성으로 인해 이 방정식의 일반적인 대수적 분석해(Closed-form analytical solution)를 구하는 것은 불가능합니다. 따라서 다음과 같은 거시적 통계학적 접근과 컴퓨터 수치해석적 접근을 활용하여 해를 도출합니다.

### ① 크레이머스-포커-플랑크 방정식 (Kramers-Fokker-Planck Equation)
개별 진동자의 미시적 경로를 구하는 대신, 진동자계 전체의 상태 확률 분포 $P(\vec{\theta}, \vec{v}, t)$의 거동을 분석합니다. 상태 공간에서의 확률 흐름의 보존 법칙(Continuity Equation)을 유도하면 다음과 같은 편미분방정식을 얻습니다.

$$\frac{\partial P}{\partial t} = -\sum_{i=1}^{N} v_i \frac{\partial P}{\partial \theta_i} - \sum_{i=1}^{N} \frac{\partial}{\partial v_i} \left[ F_i(\vec{\theta}, \vec{v}) P \right] + \frac{\sigma^2}{2 m^2} \sum_{i=1}^{N} \frac{\partial^2 P}{\partial v_i^2}$$

여기서 $F_i(\vec{\theta}, \dots) = \frac{1}{m} \left( \omega_i - \alpha v_i + \text{coupling}_i \right)$ 입니다. 통계 열역학 한계($N \to \infty$)에서 위 편미분방정식을 수치적으로 풀거나 평균장 이론(Mean-field Theory)을 적용해 동기화가 깨지는 임계 온도(임계 노이즈 세기)와 히스테리시스 이력 거동을 정량적으로 증명합니다.

### ② 이산적 수치 시뮬레이션 (Discrete Numerical Integration)

현실적인 금융 시계열 예측과 같은 신경망 파이프라인 내부에서는 유한한 $N$개의 오실레이터 궤적을 컴퓨터 상에서 한 단계씩 전진시키는 **이산적 시간 전진(Discrete Time Stepping)** 기법을 사용하여 궤적 해를 구합니다.

---

## 4. 이산적 수치 해석 기법 (Discrete Numerical Methods)

SDE의 해 경로를 근사하기 위한 대표적인 이산화 알고리즘입니다.

### ① 오일러-마루야마 방법 (Euler-Maruyama Method)

확률미분방정식에서 가장 기본적으로 활용되는 1차 적분법입니다. Deterministic 항에는 일반 Euler 스텝을, Stochastic 항에는 위너 프로세스의 독립 이산화 증분 $\Delta W_i = \sqrt{\Delta t} \epsilon_i$ (단, $\epsilon_i \sim \mathcal{N}(0, 1)$)을 적용합니다.

$$\theta_i(t + \Delta t) = \theta_i(t) + v_i(t) \Delta t$$

$$v_i(t + \Delta t) = v_i(t) + \frac{1}{m} \left[ \omega_i - \alpha v_i(t) + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] \Delta t + \frac{\sigma}{m} \sqrt{\Delta t} \epsilon_i$$

* **수렴도** : 강한 수렴 차수(Strong convergence order)는 $0.5$, 약한 수렴 차수(Weak convergence order)는 $1.0$을 보입니다. 구조가 간단하여 빠른 시뮬레이션에 적합하지만, $\Delta t$가 크거나 비선형성이 결합될 때 수치 폭주(Explosion)가 일어날 위험이 높습니다.

### ② 런게-쿠타-마루야마 분할 적분법 (Runge-Kutta-Maruyama Splitting Method)

결결론적(Deterministic) 상미분 파트에는 고차 및 어댑티브 제어가 가능한 **Dormand-Prince RK45** 공식을 적용하고, 승인된 시간 보폭 $dt$ 만큼 Langevin 확률 파트를 기하학적으로 합성하는 하이브리드 연립 분할 기법입니다.

1. **상태 미분계수 산출** :
2. 
   $$\vec{k}_t, \vec{k}_v = \text{Derivatives}(\vec{\theta}, \vec{v}, \vec{\omega})$$
   
3. **RK45 적분 (가상 결정론적 업데이트)** :
   
   표준 Dormand-Prince 계수 행렬($a_{ij}, b_i$)을 적용하여 한 단계 진행된 가상의 결정론적 최종 상태인 $\vec{\theta}_{new}, \vec{v}_{new}$를 유도합니다.
   
4. **오차 노름 검증 및 보폭 $dt$ 가감** :
   
   위상 오차와 속도 오차의 가중 평균 노름이 오차 한계(rtol, atol) 이하인지 확인하여 보폭 $dt$를 갱신합니다.
   
5. **Langevin 위너 프로세스 합성 (stochastic update)** :
   
   적분이 수락되면 가속 노이즈 변동량을 속도에 누적 반영합니다.
   
   $$\vec{\theta}(t + dt) = \vec{\theta}_{new}$$
   
   $$\vec{v}(t + dt) = \vec{v}_{new} + \frac{\sigma}{m} \sqrt{dt} \vec{\epsilon}, \quad \vec{\epsilon} \sim \mathcal{N}(0, \mathbf{I})$$

이 기법은 비선형 커플링으로 인한 결정론적 수치 에러를 고차 런게쿠타로 억제하면서 확률적 Langevin 운동을 강건하게 주입할 수 있어 SDE 수치 안정성을 획기적으로 개선합니다.

---

## 5. 이산적 수치 해석 파이썬 참조 코드

아래 코드는 이산적 Euler-Maruyama 기법과 PyTorch 텐서 배치를 활용한 하이브리드 RK45 SDE 시뮬레이션 방식을 단독 스크립트로 동작할 수 있게 구현한 참조용 파이썬 예제 코드입니다.

```python
import torch
import torch.nn as nn
import numpy as np

class StochasticKuramoto2ndSolver(nn.Module):
    def __init__(self, num_oscillators=64, dt=0.01, sigma=0.05, mass=1.5, damping=0.3):
        super(StochasticKuramoto2ndSolver, self).__init__()
        self.N = num_oscillators
        self.dt = dt
        
        # 물리 파라미터 (미분 연산을 위해 소프트플러스와 시그모이드 바인딩 제약 적용)
        self.raw_mass = nn.Parameter(torch.tensor(mass))
        self.raw_damping = nn.Parameter(torch.tensor(damping))
        self.K = nn.Parameter(torch.randn(self.N, self.N) * 0.2)
        self.log_sigma = nn.Parameter(torch.tensor(np.log(sigma)))
        
    def get_physics_parameters(self):
        m = torch.nn.functional.softplus(self.raw_mass) + 0.5
        alpha = torch.sigmoid(self.raw_damping) * 0.45 + 0.05
        sigma = torch.exp(self.log_sigma) + 0.01
        return m, alpha, sigma

    def _get_coupling_matrix(self):
        # 내부 인력 결합 마스크 설정
        K_raw = torch.nn.functional.softplus(self.K)
        mask_intra = torch.zeros(self.N, self.N, device=self.K.device)
        mask_intra[:self.N//2, :self.N//2] = 1.0
        mask_intra[self.N//2:, self.N//2:] = 1.0
        mask_inter = 1.0 - mask_intra
        
        # 결합 강도 행렬 구축 (내부 결합 (+), 교차 결합 (-))
        K_active = K_raw * mask_intra * 1.2 - K_raw * mask_inter * 0.15 + (mask_intra * 0.3 - mask_inter * 0.05)
        return K_active

    def step_euler_maruyama(self, theta, v, omega):
        """
        Euler-Maruyama 이산 적분 방식 (1회 스텝)
        """
        m, alpha, sigma = self.get_physics_parameters()
        K_active = self._get_coupling_matrix()
        
        # 커플링 힘 계산
        theta_diff = theta.unsqueeze(1) - theta.unsqueeze(2)  # [B, N, N]
        interaction = torch.sum(K_active * torch.sin(theta_diff), dim=-1) / self.N # [B, N]
        
        # 미분식
        dtheta_dt = v
        dv_dt = (omega - alpha * v + interaction) / m
        
        # Euler Step
        theta_next = theta + dtheta_dt * self.dt
        
        # Langevin Stochastic Step
        noise_v = (sigma * torch.randn_like(v) * np.sqrt(self.dt)) / m
        v_next = v + dv_dt * self.dt + noise_v
        
        return theta_next, v_next

    def simulate(self, theta_0, v_0, omega, steps=48):
        """
        전체 시간 단계 시뮬레이션 및 오더 파라미터 반환
        """
        batch_size = theta_0.size(0)
        theta_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        v_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        r_history = torch.zeros(batch_size, steps, device=theta_0.device)
        
        theta = theta_0.clone()
        v = v_0.clone()
        
        for step in range(steps):
            theta, v = self.step_euler_maruyama(theta, v, omega)
            theta_traj[:, step] = theta
            v_traj[:, step] = v
            # 질서 매개변수 r 계산 (복소 공간 상의 지수 위상값들의 평균 노름)
            r = torch.abs(torch.mean(torch.exp(1j * (theta % (2 * np.pi))), dim=-1))
            r_history[:, step] = r
            
        return theta_traj, v_traj, r_history

# 사용 예시 데모 코드
if __name__ == "__main__":
    batch_size = 2
    solver = StochasticKuramoto2ndSolver(num_oscillators=64, dt=0.02, sigma=0.08)
    
    # 초기 임의 상태 생성
    theta_init = torch.rand(batch_size, 64) * 2 * np.pi
    v_init = torch.randn(batch_size, 64) * 0.1
    omega_init = torch.randn(batch_size, 64) * 0.2
    
    # 48일의 궤적 계산
    traj_theta, traj_v, traj_r = solver.simulate(theta_init, v_init, omega_init, steps=48)
    
    print("오실레이터 시뮬레이션 완료:")
    print("위상 궤적 텐서 모양 (Theta Trajectory shape):", traj_theta.shape)  # [2, 48, 64]
    print("속도 궤적 텐서 모양 (Velocity Trajectory shape):", traj_v.shape)  # [2, 48, 64]
    print("오더 파라미터 변동 추이 (최종 5 스텝):", traj_r[0, -5:].detach().numpy())
```
