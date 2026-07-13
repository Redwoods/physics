# 결정론적 1차 쿠라모토(Kuramoto) 모형의 이론 및 수치적 해석 기법

본 보고서는 관성(질량)과 확률적 외부 잡음이 고려되지 않은 순수 결정론적 비선형 시스템인 **'결정론적 1차 쿠라모토 모형(Deterministic 1st-order Kuramoto Model)'**의 지배 방정식, 수학적 해석 방법 및 컴퓨터로 이산적인 해를 구하는 수치 해석 알고리즘과 참조 파이썬 코드를 상세히 정리한 문서입니다.

---

## 1. 결정론적 1차 쿠라모토 지배 방정식

결정론적 1차 쿠라모토 모델은 마찰 저항이 극도로 큰 과감쇠 한계(Overdamped limit) 상태를 기술하며, 오실레이터들의 위상 변화율(회전 속도)이 위상차에 따른 결합력에 의해 지연 시간 없이 즉각 결정된다고 가정합니다. 지배 상미분방정식(ODE)은 다음과 같이 표현됩니다.

$$\frac{d\theta_i(t)}{dt} = \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t))$$

이 식에서 각 물리량의 정의는 다음과 같습니다:
* $\theta_i(t)$: 오실레이터 $i$의 위상 (Phase)
* $\omega_i$: 오실레이터 $i$의 고유 진동수 (Intrinsic Frequency)
* $K$: 결합 강도 (Coupling strength)
* $A_{ij}$: 인접 행렬(Adjacency Matrix, 네트워크 연결 가중치 및 마스크)
* $N$: 전체 오실레이터의 수

---

## 2. 평균장 단순화와 자가일치 조건 (Mean-field Formulation)

결정론적 쿠라모토 모델 역시 집단적인 동기화 정도를 측정하기 위해 복소 **질서 매개변수(Order Parameter, $r$)**와 집단 평균 위상($\psi$)을 다음과 같이 정의합니다.

$$r(t) e^{i \psi(t)} = \frac{1}{N} \sum_{j=1}^{N} e^{i \theta_j(t)}$$

이 관계를 복소수 평면 상에서 전개하면 $\frac{1}{N} \sum_j \sin(\theta_j - \theta_i) = r \sin(\psi - \theta_i)$ 임을 쉽게 증명할 수 있습니다. 이를 지배 방정식에 대입하면 각각의 오실레이터가 나머지 $N-1$개의 진동자들과 일일이 연산하는 대신, 공통의 평균장 $r(t)$에만 동적으로 반응하는 평균장 방정식으로 단순화됩니다:

$$\frac{d\theta_i(t)}{dt} = \omega_i + K r(t) \sin(\psi(t) - \theta_i(t))$$

---

## 3. 수학적 해를 구하는 방법론 (Mathematical Solution Methods)

오실레이터 수가 무한대($N \to \infty$)이고 시스템이 정상 동기화 상태(Stationary synchronized state)에 도달했다고 가정합니다. 이때 평균 위상 $\psi(t)$는 일정한 집단 평균 진동수 $\Omega$로 회전하며 ($\psi(t) = \Omega t$), 질서 매개변수 $r$은 시간에 따라 일정한 상수 값을 유지합니다.

평균 위상과 함께 회전하는 회전 좌표계를 도입하여 상대 위상차를 $\phi_i = \theta_i - \Omega t$로 정의하면, 방정식은 자율 상미분방정식(Autonomous ODE) 형태로 바뀝니다:

$$\frac{d\phi_i}{dt} = \omega_i - \Omega - K r \sin\phi_i$$

이 좌표계에서 오실레이터들의 거동은 고유 진동수 크기에 따라 두 부류로 명확히 분리됩니다.

### ① 동기화 진동자 (Phase-locked Oscillators)
고유 진동수가 평균 진동수 근방에 위치하여 결합력 장벽보다 작거나 같은 진동자들입니다.
$$|\omega_i - \Omega| \le K r$$
이들은 회전 좌표계 상에서 특정 고정 위상각(Fixed point)에 수렴하여 상대 속도가 0이 됩니다 ($\frac{d\phi_i}{dt} = 0$). 위상 평형 조건은 다음과 같습니다.
$$\sin\phi_i = \frac{\omega_i - \Omega}{K r}$$

### ② 비동기화 표류 진동자 (Drifting Oscillators)
고유 진동수가 결합력 한계를 초과하여 집단 속도에 결속되지 못하고 독립적으로 계속 표류하는 진동자들입니다.
$$|\omega_i - \Omega| > K r$$
이들은 동기화 상태에 도달하지 못하고 계속 360도 회전하지만, 결합력의 영향으로 인해 특정 위상각 부근을 지날 때 회전 속도가 느려지고 반대쪽에서는 빨라집니다. 정상 상태에서 이 표류 오실레이터들이 임의의 각도 $\phi$에 존재할 정밀 확률 밀도 $n(\phi, \omega)$는 속도 크기에 반비례합니다:
$$n(\phi, \omega) = \frac{C}{|d\phi/dt|} = \frac{\sqrt{(\omega - \Omega)^2 - (Kr)^2}}{2\pi |\omega - \Omega - Kr \sin\phi|}$$

### ③ 쿠라모토 자가일치 관계식 및 임계값 ($K_c$)
동기화된 진동자들과 표류하는 진동자들의 기여도를 합산하여 복소 질서 매개변수의 자가일치 조건(Self-consistency relation) 적분식을 세우면 다음과 같습니다. 고유 진동수 분포 $g(\omega)$가 우함수(symmetric)이고 단봉(unimodal) 분포일 때 집단 진동수 $\Omega$는 분포의 중심과 일치하고, 질서 매개변수 $r$은 다음 비선형 방정식을 만족합니다.

$$r = K r \int_{-\pi/2}^{\pi/2} \cos^2\phi \, g(Kr \sin\phi) d\phi$$

질서 매개변수가 비동기화 상태($r=0$)에서 동기화 상태($r > 0$)로 상전이를 일으키는 임계 결합 강도 $K_c$는 $r \to 0^+$ 극한을 취하여 적분함으로써 수학적으로 정확히 유도됩니다.

$$K_c = \frac{2}{\pi g(0)}$$

예를 들어, 고유 진동수 분포 $g(\omega)$가 표준 편차 $\gamma$를 가지는 Lorentzian 분포인 경우 임계값은 $K_c = 2\gamma$가 되며, $K > K_c$ 일 때 생성되는 정상 상태 질서 매개변수 $r$의 분기 거동은 다음과 같은 멱법칙을 따릅니다:
$$r \approx \sqrt{\frac{16}{\pi K_c^3} \frac{K - K_c}{-g''(0)}} \propto (K - K_c)^{1/2}$$

---

## 4. 이산적 수치 해석 기법 (Discrete Numerical Methods)

컴퓨터 상에서 1차 결정론적 ODE의 궤적을 전진시키는 대표적인 수치 적분 해법입니다.

### ① Euler 방법 (1차 전진 오일러)
가장 단순한 시간 전진 적분 기법입니다.
$$\theta_i(t + \Delta t) = \theta_i(t) + \left[ \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] \Delta t$$
연산 속도가 극단적으로 빠르지만 수치 안정성을 담보하기 위해 보폭 $\Delta t$를 매우 작게 설정해야 합니다.

### ② Runge-Kutta 4차 방법 (RK4)
상미분방정식에서 보편적으로 사용되는 높은 정확도의 이산화 기법입니다. 보폭 $\Delta t$를 크게 잡더라도 4차 오차 정확도($\mathcal{O}(\Delta t^4)$)로 안정적으로 해 궤적을 계산합니다.
$$k_1 = f(\theta(t), t)$$
$$k_2 = f\left(\theta(t) + \frac{\Delta t}{2} k_1, t + \frac{\Delta t}{2}\right)$$
$$k_3 = f\left(\theta(t) + \frac{\Delta t}{2} k_2, t + \frac{\Delta t}{2}\right)$$
$$k_4 = f(\theta(t) + \Delta t k_3, t + \Delta t)$$
$$\theta(t + \Delta t) = \theta(t) + \frac{\Delta t}{6} (k_1 + 2k_2 + 2k_3 + k_4)$$

---

## 5. 이산적 수치 해석 파이썬 참조 코드

아래 코드는 4차 Runge-Kutta (RK4) 기법과 PyTorch 텐서 배치 연산을 결합하여 1차 결정론적 쿠라모토 모형의 거동을 시뮬레이션하는 파이썬 코드 예제입니다.

```python
import torch
import torch.nn as nn

class DeterministicKuramoto1stSolver(nn.Module):
    def __init__(self, num_oscillators=64, dt=0.02):
        super(DeterministicKuramoto1stSolver, self).__init__()
        self.N = num_oscillators
        self.dt = dt
        
        # 결합 강도 매트릭스 파라미터화 (softplus 바인딩 적용)
        self.K = nn.Parameter(torch.randn(self.N, self.N) * 0.2)
        
    def _derivatives(self, theta, omega):
        """1계 미분계수 산출: dtheta_dt = omega + interaction"""
        K_raw = torch.nn.functional.softplus(self.K)
        
        # 다중 물리 블록 결합 강도 설정 (내부 인력, Inter-group 척력)
        mask_intra = torch.zeros(self.N, self.N, device=self.K.device)
        mask_intra[:self.N//2, :self.N//2] = 1.0
        mask_intra[self.N//2:, self.N//2:] = 1.0
        mask_inter = 1.0 - mask_intra
        K_active = K_raw * mask_intra * 1.2 - K_raw * mask_inter * 0.15 + (mask_intra * 0.3 - mask_inter * 0.05)
        
        theta_diff = theta.unsqueeze(1) - theta.unsqueeze(2)  # [B, N, N]
        interaction = torch.sum(K_active * torch.sin(theta_diff), dim=-1) / self.N # [B, N]
        
        dtheta_dt = omega + interaction
        return dtheta_dt

    def step_rk4(self, theta, omega):
        """
        4차 Runge-Kutta (RK4) 1회 스텝
        """
        k1 = self._derivatives(theta, omega)
        
        theta_2 = theta + 0.5 * self.dt * k1
        k2 = self._derivatives(theta_2, omega)
        
        theta_3 = theta + 0.5 * self.dt * k2
        k3 = self._derivatives(theta_3, omega)
        
        theta_4 = theta + self.dt * k3
        k4 = self._derivatives(theta_4, omega)
        
        theta_next = theta + (self.dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        return theta_next

    def simulate(self, theta_0, omega, steps=48):
        """
        결정론적 다단계 시간 적분 및 오더 파라미터 기록
        """
        batch_size = theta_0.size(0)
        theta_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        r_history = torch.zeros(batch_size, steps, device=theta_0.device)
        
        theta = theta_0.clone()
        
        for step in range(steps):
            theta = self.step_rk4(theta, omega)
            theta_traj[:, step] = theta
            
            # 오더 파라미터 r 계산
            r = torch.abs(torch.mean(torch.exp(1j * (theta % (2 * np.pi))), dim=-1))
            r_history[:, step] = r
            
        return theta_traj, r_history

# 사용 예시 데모 스크립트
if __name__ == "__main__":
    batch_size = 2
    solver = DeterministicKuramoto1stSolver(num_oscillators=64, dt=0.02)
    
    # 초기 위상 및 고유 진동수 무작위 설정
    theta_init = torch.rand(batch_size, 64) * 2 * np.pi
    omega_init = torch.randn(batch_size, 64) * 0.15
    
    # 48일 시뮬레이션 전개
    traj_theta, traj_r = solver.simulate(theta_init, omega_init, steps=48)
    
    print("결정론적 1차 오실레이터 RK4 시뮬레이션 완료:")
    print("위상 궤적 텐서 모양 (Theta Trajectory shape):", traj_theta.shape)  # [2, 48, 64]
    print("최종 5일의 오더 파라미터 r 값:", traj_r[0, -5:].detach().numpy())
```
