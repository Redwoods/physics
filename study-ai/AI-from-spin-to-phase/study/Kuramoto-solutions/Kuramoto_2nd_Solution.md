# 결정론적 2차 쿠라모토(Kuramoto) 모형의 이론 및 수치적 해석 기법

본 보고서는 질량(관성)과 감쇠(저항)가 존재하며 확률적 노이즈가 배제된 결정론적 상미분 동역학계 모델인 **'결정론적 2차 쿠라모토 모형(Deterministic 2nd-order Kuramoto Model)'**의 지배 방정식, 수학적 해석 방법 및 컴퓨터로 이산적인 해를 구하는 수치 해석 알고리즘과 참조 파이썬 코드를 상세히 정리한 문서입니다.

---

## 1. 결정론적 2차 쿠라모토 지배 방정식

2차 쿠라모토 모델은 개별 오실레이터에 물리적인 회전 **질량(관성, $m$)**과 저항 역할을 하는 **감쇠 계수(마찰, $\alpha$)**를 추가하여 뉴턴 역학적인 가속도 운동을 기술합니다. 지배 방정식은 다음과 같은 2계 상미분방정식(2nd-order ODE)으로 정의됩니다.

$$m \ddot{\theta}_i(t) + \alpha \dot{\theta}_i(t) = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t))$$

이 식에서 각 물리량의 정의는 다음과 같습니다:
* $\theta_i(t)$: 오실레이터 $i$의 위상 (Phase)
* $v_i(t) = \dot{\theta}_i(t)$: 오실레이터 $i$의 각속도 (Angular Velocity)
* $m$: 오실레이터의 질량 (관성 크기, Inertia)
* $\alpha$: 감쇠 계수 (마찰/저항 계수, Damping coefficient)
* $\omega_i$: 오실레이터 $i$의 고유 진동수 (Intrinsic Frequency)
* $K$: 결합 강도 (Coupling strength)
* $A_{ij}$: 인접 행렬(Adjacency Matrix, 마스크 가중치)
* $N$: 전체 오실레이터 수

---

## 2. 상태 공간 표현식 및 평균장 단순화 (State-Space & Mean-Field)

컴퓨터 상에서 수치 연산을 효과적으로 수행하기 위해 속도 상태 변수 $v_i(t) = \dot{\theta}_i(t)$를 도입하여 2계 ODE를 연립 1계 ODE 시스템으로 분리합니다.

$$\frac{d\theta_i(t)}{dt} = v_i(t)$$

$$\frac{dv_i(t)}{dt} = \frac{1}{m} \left[ \omega_i - \alpha v_i(t) + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right]$$

또한, 집단의 평균 흐름을 대변하는 복소 **질서 매개변수(Order Parameter, $r$)**와 집단 평균 위상($\psi$)을 다음과 같이 정의합니다:

$$r(t) e^{i \psi(t)} = \frac{1}{N} \sum_{j=1}^{N} e^{i \theta_j(t)}$$

이를 결합하면 각 오실레이터의 질량 가속 운동은 공유하는 평균장 $r(t)$에만 동기화 힘을 받는 형태로 단순화됩니다:

$$m \ddot{\theta}_i(t) + \alpha \dot{\theta}_i(t) = \omega_i + K r(t) \sin(\psi(t) - \theta_i(t))$$

---

## 3. 수학적 해를 구하는 방법론 (Mathematical Solution Methods)

오실레이터의 개수가 충분히 많고 시스템이 정상 상태에 도달해 평균 위상 $\psi(t)$가 일정한 평균 진동수 $\Omega$로 회전하며 질서 매개변수 $r$이 일정한 상수로 안착했다고 가정합니다.

평균 위상의 속도로 회전하는 동기식 회전 좌표계를 적용하여 상대 위상차를 $\phi_i = \theta_i - \Omega t$로 두면 다음과 같은 지배 방정식이 유도됩니다:

$$m \ddot{\phi}_i + \alpha \dot{\phi}_i = \omega_i - \alpha \Omega - K r \sin\phi_i$$

### ① 동기화 진동자 (Phase-locked Oscillators)
상대 속도와 상대 가속도가 완전히 0이 되어 특정 평형 위상각에 고정되는 오실레이터들입니다 ($\dot{\phi}_i = \ddot{\phi}_i = 0$). 위상 고정 평형 조건식은 다음과 같습니다.

$$\sin\phi_i = \frac{\omega_i - \alpha \Omega}{K r}$$

이 평형 실근이 존재하기 위한 필요충분조건은 고유진동수의 편차가 결합력보다 작아야 함을 의미합니다:

$$|\omega_i - \alpha \Omega| \le K r$$

이때 1차 모델과 달리 2차 모델은 동일 범위 내에 두 개의 평형 상태(안정 노드와 안장점)가 공존하며, $\cos\phi_i > 0$ 을 만족하는 위상 각도만 국소적 안정성(Local stability)을 확보합니다.

### ② 비동기화 진동자 및 표류 (Drifting Oscillators)
고유진동수 편차가 결합 한계를 벗어나 동기화에 잡히지 않고 계속 회전하는 진동자들입니다 ($|\omega_i - \alpha \Omega| > K r$). 이들은 마찰력과 질량 관성의 이중 제어를 받으며 끊임없이 주파수가 요동치는 표류 운동을 전개합니다.

### ③ 1차 불연속 상전이 및 히스테리시스 이력 현상 (Hysteresis)
2차 쿠라모토 모델의 가장 핵심적인 수학적·물리적 특성은 1차 모델의 연속 상전이와 달리 **1차 불연속 상전이(1st-order Discontinuous Phase Transition)**와 **히스테리시스(이력 현상)**를 유발한다는 점입니다.
* **불연속 상전이**: 결합 강도 $K$를 점진적으로 늘릴 때, 동기화 지표 $r$이 0에서 완만하게 상승하지 않고 임계점 $K_c$에 도달하는 순간 갑자기 유한한 동기화 값으로 **불연속적인 점프(Discontinuous Jump)**를 일으킵니다.

* **히스테리시스 루프**: 이미 동기화된 시스템($r > 0$)의 결합도를 다시 약화시킬 때, 결합도가 $K_c$보다 작아지더라도 동기화가 즉시 파괴되지 않고 더 낮은 하한 임계점 $K_c^-$까지 동기화 상태가 잔존합니다. 즉, $K \in [K_c^-, K_c]$ 구간에서 시스템은 과거의 경로에 따라 동기화 상태와 비동기화 상태 중 하나로 결정되는 다중 안정성(Multistability) 및 시간적 이력을 가집니다. 이 관성(질량)에 의한 히스테리시스 특성은 주식 시장에서 정보 반영이 지연되며 급등락 추세가 지속되는 오버슈팅 현상을 설명하는 물리적 근거가 됩니다.

---

## 4. 이산적 수치 해석 기법 (Discrete Numerical Methods)

2차 상미분방정식을 이산적으로 해결하기 위해 상태 공간 텐서 $\vec{x}(t) = [\vec{\theta}(t), \vec{v}(t)]^T$를 구성하고 시간 스텝 $\Delta t$ 단위로 전진시킵니다.

### ① Euler 방법 (1차 연립 적분)
가장 원초적인 이산화 스텝 기법입니다.

$$\theta_i(t + \Delta t) = \theta_i(t) + v_i(t) \Delta t$$

$$v_i(t + \Delta t) = v_i(t) + \frac{1}{m} \left[ \omega_i - \alpha v_i(t) + \text{coupling}_i(t) \right] \Delta t$$

### ② Runge-Kutta 4차 방법 (RK4)
상태 공간 연립 미분 벡터 방정식 $\vec{x}' = \vec{f}(\vec{x}, t)$에 대해 오차 차수 $\mathcal{O}(\Delta t^4)$를 제공하는 표준 수치해석 기법입니다.

$$\vec{k}_1 = \vec{f}(\vec{x}(t), t)$$

$$\vec{k}_2 = \vec{f}\left(\vec{x}(t) + \frac{\Delta t}{2} \vec{k}_1, t + \frac{\Delta t}{2}\right)$$

$$\vec{k}_3 = \vec{f}\left(\vec{x}(t) + \frac{\Delta t}{2} \vec{k}_2, t + \frac{\Delta t}{2}\right)$$

$$\vec{k}_4 = \vec{f}(\vec{x}(t) + \Delta t \vec{k}_3, t + \Delta t)$$

$$\vec{x}(t + \Delta t) = \vec{x}(t) + \frac{\Delta t}{6} (\vec{k}_1 + 2\vec{k}_2 + 2\vec{k}_3 + \vec{k}_4)$$


---

## 5. 이산적 수치 해석 파이썬 참조 코드

아래 코드는 PyTorch 배치 연산과 4차 Runge-Kutta (RK4) 수치 적분을 이용하여 결정론적 2차 쿠라모토 모형의 궤적을 산출하는 단독 구동용 파이썬 클래스 예제입니다.

```python
import torch
import torch.nn as nn

class DeterministicKuramoto2ndSolver(nn.Module):
    def __init__(self, num_oscillators=64, dt=0.02, mass=1.5, damping=0.3):
        super(DeterministicKuramoto2ndSolver, self).__init__()
        self.N = num_oscillators
        self.dt = dt
        
        # 학습 및 최적화 가능한 물리 파라미터 제약화 설정
        self.raw_mass = nn.Parameter(torch.tensor(mass))
        self.raw_damping = nn.Parameter(torch.tensor(damping))
        self.K = nn.Parameter(torch.randn(self.N, self.N) * 0.2)
        
    def get_physics_parameters(self):
        m = torch.nn.functional.softplus(self.raw_mass) + 0.5
        alpha = torch.sigmoid(self.raw_damping) * 0.45 + 0.05
        return m, alpha

    def _derivatives(self, theta, v, omega):
        """2차 연립 미분계수 산출: dtheta_dt = v, dv_dt = (omega - alpha*v + interaction)/m"""
        m, alpha = self.get_physics_parameters()
        K_raw = torch.nn.functional.softplus(self.K)
        
        # 다중 물리 내부/그룹간 결합 설정
        mask_intra = torch.zeros(self.N, self.N, device=self.K.device)
        mask_intra[:self.N//2, :self.N//2] = 1.0
        mask_intra[self.N//2:, self.N//2:] = 1.0
        mask_inter = 1.0 - mask_intra
        K_active = K_raw * mask_intra * 1.2 - K_raw * mask_inter * 0.15 + (mask_intra * 0.3 - mask_inter * 0.05)
        
        theta_diff = theta.unsqueeze(1) - theta.unsqueeze(2)  # [B, N, N]
        interaction = torch.sum(K_active * torch.sin(theta_diff), dim=-1) / self.N # [B, N]
        
        # dtheta_dt 및 dv_dt 반환
        dtheta_dt = v
        dv_dt = (omega - alpha * v + interaction) / m
        return dtheta_dt, dv_dt

    def step_rk4(self, theta, v, omega):
        """
        위상과 속도의 연립 4차 Runge-Kutta (RK4) 이산화 스텝
        """
        k1_t, k1_v = self._derivatives(theta, v, omega)
        
        t2 = theta + 0.5 * self.dt * k1_t
        v2 = v + 0.5 * self.dt * k1_v
        k2_t, k2_v = self._derivatives(t2, v2, omega)
        
        t3 = theta + 0.5 * self.dt * k2_t
        v3 = v + 0.5 * self.dt * k2_v
        k3_t, k3_v = self._derivatives(t3, v3, omega)
        
        t4 = theta + self.dt * k3_t
        v4 = v + self.dt * k3_v
        k4_t, k4_v = self._derivatives(t4, v4, omega)
        
        theta_next = theta + (self.dt / 6.0) * (k1_t + 2*k2_t + 2*k3_t + k4_t)
        v_next     = v + (self.dt / 6.0) * (k1_v + 2*k2_v + 2*k3_v + k4_v)
        
        return theta_next, v_next

    def simulate(self, theta_0, v_0, omega, steps=48):
        """
        결정론적 다단계 시간 연립 시뮬레이션 수행 및 질서 매개변수 r 출력
        """
        batch_size = theta_0.size(0)
        theta_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        v_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        r_history = torch.zeros(batch_size, steps, device=theta_0.device)
        
        theta = theta_0.clone()
        v = v_0.clone()
        
        for step in range(steps):
            theta, v = self.step_rk4(theta, v, omega)
            theta_traj[:, step] = theta
            v_traj[:, step] = v
            
            # 오더 파라미터 r 계산
            r = torch.abs(torch.mean(torch.exp(1j * (theta % (2 * np.pi))), dim=-1))
            r_history[:, step] = r
            
        return theta_traj, v_traj, r_history

# 데모 실행부
if __name__ == "__main__":
    batch_size = 2
    solver = DeterministicKuramoto2ndSolver(num_oscillators=64, dt=0.02, mass=1.8, damping=0.25)
    
    # 초기 상태 및 속도 설정
    theta_init = torch.rand(batch_size, 64) * 2 * np.pi
    v_init = torch.randn(batch_size, 64) * 0.05
    omega_init = torch.randn(batch_size, 64) * 0.1
    
    # 48 스텝 시뮬레이션 전개
    traj_theta, traj_v, traj_r = solver.simulate(theta_init, v_init, omega_init, steps=48)
    
    m_val, alpha_val = solver.get_physics_parameters()
    print(f"물리 파라미터 확인 -> 질량: {m_val.item():.4f}, 감쇠비: {alpha_val.item():.4f}")
    print("결정론적 2차 오실레이터 연립 RK4 시뮬레이션 완료:")
    print("위상 궤적 텐서 모양 (Theta Trajectory shape):", traj_theta.shape)  # [2, 48, 64]
    print("속도 궤적 텐서 모양 (Velocity Trajectory shape):", traj_v.shape)  # [2, 48, 64]
    print("최종 5일의 오더 파라미터 r 값:", traj_r[0, -5:].detach().numpy())
```
