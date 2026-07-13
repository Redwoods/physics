# 확률적 1차 쿠라모토(Kuramoto) SDE 모형의 이론 및 수치적 해석 기법

본 보고서는 관성(질량)이 없는 1차 쿠라모토 진동자계에 확률적 노이즈가 도입된 **'확률적 1차 쿠라모토 SDE 모형(Stochastic 1st-order Kuramoto SDE Model)'**의 지배 방정식, 수학적 해석 방법 및 컴퓨터로 이산적인 해를 구하는 수치 해석 알고리즘과 참조 파이썬 코드를 상세히 정리한 문서입니다.

---

## 1. 확률적 1차 쿠라모토 SDE 지배 방정식

1차 쿠라모토 모델은 오실레이터가 상호작용할 때 지연 시간 없이 위상을 즉각적으로 조정한다고 가정합니다. 오실레이터들의 위상 상태 갱신에 외부 노이즈(시장 잡음 및 확률적 교란)를 나타내는 **표준 위너 프로세스(Wiener Process)**를 가산 형태로 결합하면 다음과 같은 이토 확률미분방정식(Itô SDE)을 얻습니다.

$$d\theta_i(t) = \left[ \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] dt + \sigma dW_i(t)$$

이 식에서 각 물리량의 정의는 다음과 같습니다:
* $\theta_i(t)$: 오실레이터 $i$의 위상 (Phase)
* $\omega_i$: 오실레이터 $i$의 고유 진동수 (Intrinsic Frequency)
* $K$: 결합 강도 (Coupling strength)
* $A_{ij}$: 인접 행렬(Adjacency Matrix, 마스크 가중치)
* $\sigma$: 노이즈 강도 (Noise intensity)
* $W_i(t)$: 표준 위너 프로세스 (Wiener Process, Brownian Motion)로, 이산화 증분은 $dW_i(t) \approx \mathcal{N}(0, dt)$를 따릅니다.

상태 변수인 위상 $\theta_i$에 잡음이 직접 더해지는 **덧셈 노이즈(Additive Noise)** 형태이므로, 이토(Itô) 적분과 스트라토노비치(Stratonovich) 적분의 수학적 형태가 완벽히 일치하여 수식 해석이 직관적입니다.

---

## 2. 평균장 이론과 단순화 (Mean-field Approximation)

1차 쿠라모토 모델은 시스템 전체의 동기화 정도를 정량화하기 위해 **질서 매개변수(Order Parameter, $r$)**와 집단의 평균 위상($\psi$)을 다음과 같이 정의합니다.

$$r(t) e^{i \psi(t)} = \frac{1}{N} \sum_{j=1}^{N} e^{i \theta_j(t)}$$

여기서 질서 매개변수 $r(t) \in [0, 1]$는 0에 가까우면 모든 오실레이터가 무질서한 비동기화 상태임을 의미하고, 1에 가까우면 완벽한 동기화 상태에 도달했음을 의미합니다.

위의 질서 매개변수의 정의를 원래의 SDE 지배 방정식에 대입하여 삼각함수 공식을 정리하면, $N^2$번의 상호작용 계산을 거치지 않고 각 오실레이터가 평균장에만 반응하는 형태의 **평균장 지배 방정식(Mean-field Langevin Equation)**으로 단순화할 수 있습니다:

$$d\theta_i(t) = \left[ \omega_i + K r(t) \sin(\psi(t) - \theta_i(t)) \right] dt + \sigma dW_i(t)$$

이 식은 복잡한 다체 상호작용 문제를 평균장 $r(t)$ 하의 1체 운동으로 치환할 수 있어 수학적 해석 및 수치 연산 효율성을 극대화합니다.

---

## 3. 수학적 해를 구하는 방법론 (Mathematical Solution Methods)

사인 결합 항의 비선형성으로 인해 개별 궤적의 일반적인 analytical solution을 구하는 것은 불가능하지만, 통계역학적으로 다음과 같은 해석 기법을 사용합니다.

### ① 쿠라모토-사카구치 포커-플랑크 방정식 (Kuramoto-Sakaguchi Fokker-Planck Equation)
오실레이터의 개수가 무한대($N \to \infty$)로 수렴하는 연속체 한계(Thermodynamic limit)에서, 특정 고유진동수 $\omega$를 지닌 오실레이터가 시점 $t$에 위상 $\theta$에 존재할 확률 밀도 함수 $P(\theta, \omega, t)$의 거동은 다음 편미분방정식을 따릅니다.

$$\frac{\partial P}{\partial t} = -\frac{\partial}{\partial \theta} \left[ \left( \omega + K r \sin(\psi - \theta) \right) P \right] + \frac{\sigma^2}{2} \frac{\partial^2 P}{\partial \theta^2}$$

이 방정식에서 정상 상태($\partial P / \partial t = 0$) 조건과 질서 매개변수 $r$의 자가일치 조건(Self-consistency relation)을 유도하면 다음과 같이 동기화 상전이(Phase Transition)가 일어나는 임계 결합 강도 $K_c$를 노이즈 강도 $\sigma$의 함수로 증명할 수 있습니다. 고유 진동수 분포 $g(\omega)$가 로렌츠 분포(Lorentzian Distribution)이고 폭이 $\gamma$일 때:

$$K_c = 2 \left( \gamma + \frac{\sigma^2}{2} \right)$$

즉, 시스템에 존재하는 확률적 잡음 강도 $\sigma$가 클수록 오실레이터들을 정렬하기 위해 필요한 최소 결합 강도 $K_c$의 문턱(Threshold)이 선형적으로 높아짐을 수학적으로 알 수 있습니다.

---

## 4. 이산적 수치 해석 기법 (Discrete Numerical Methods)

컴퓨터 연산을 위해 시계열 시간 축을 이산적인 시간 보폭 $\Delta t$ 단위로 나누어 궤적을 전진시키는 해법들입니다.

### ① 오일러-마루야마 방법 (Euler-Maruyama Method)
1차 SDE 모형에 Euler-Maruyama를 적용하면 가속도가 없는 심플한 형태로 이산화됩니다.

$$\theta_i(t + \Delta t) = \theta_i(t) + \left[ \omega_i + \frac{K}{N} \sum_{j=1}^{N} A_{ij} \sin(\theta_j(t) - \theta_i(t)) \right] \Delta t + \sigma \sqrt{\Delta t} \epsilon_i$$

여기서 $\epsilon_i \sim \mathcal{N}(0, 1)$은 가우시안 정규분포 난수입니다. 구현이 단순하여 연산 효율이 매우 뛰어납니다.

### ② 런게-쿠타-마루야마 분할 적분법 (Runge-Kutta-Maruyama Splitting Method)
결정론적 동역학 항은 고차 오차 정밀도와 어댑티브 보폭을 제공하는 **Dormand-Prince RK45** 적분기로 가상의 $\theta_{new}$를 유도하고, 수락된 시간 보폭 $dt$에 대해 Langevin 잡음을 더해주는 분할 결합 기법입니다. [hkino_Kuramoto_1st_SDE.py](file:///home/redwoods/gdrive/AI/AI_ALL/New_methods/Phase_Deep_Learning/Deep_Research/Z_Stock_Prediction/hkino_Kuramoto_1st_SDE.py)에 실제 채택된 수치 적분 구조입니다.

1. **결정론적 궤적 계산**:
   $$\theta_{new} = \text{RK45\_Deterministic\_Update}(\theta, \omega, dt)$$
2. **Wiener Step 가산**:
   수락 마스크가 활성화되면 오차 요건에 맞추어 보정된 $dt$ 크기의 난수 잡음을 갱신된 위상에 가산합니다.
   $$\theta(t + dt) = \theta_{new} + \sigma \sqrt{dt} \epsilon, \quad \epsilon \sim \mathcal{N}(0, \mathbf{I})$$

이 방법은 비선형 결합력의 시간 누적 오차를 고차 RK45로 모니터링하여 보폭을 유연하게 잡고, 확률적 잡음을 결합하므로 안정적이고 정확한 시뮬레이션 경로를 제공합니다.

---

## 5. 이산적 수치 해석 파이썬 참조 코드

아래 코드는 Euler-Maruyama 기법과 PyTorch 배치 텐서 연산을 이용하여 1차 확률적 쿠라모토 SDE를 해결하는 참조용 파이썬 예제 코드입니다.

```python
import torch
import torch.nn as nn
import numpy as np

class StochasticKuramoto1stSolver(nn.Module):
    def __init__(self, num_oscillators=64, dt=0.01, sigma=0.05):
        super(StochasticKuramoto1stSolver, self).__init__()
        self.N = num_oscillators
        self.dt = dt
        
        # 결합 파라미터 및 학습 가능한 노이즈 세기 파라미터 (exp(log_sigma)로 양수 보장)
        self.K = nn.Parameter(torch.randn(self.N, self.N) * 0.2)
        self.log_sigma = nn.Parameter(torch.tensor(np.log(sigma)))
        
    def get_physics_parameters(self):
        sigma = torch.exp(self.log_sigma) + 0.01
        return sigma

    def _get_coupling_matrix(self):
        K_raw = torch.nn.functional.softplus(self.K)
        mask_intra = torch.zeros(self.N, self.N, device=self.K.device)
        mask_intra[:self.N//2, :self.N//2] = 1.0
        mask_intra[self.N//2:, self.N//2:] = 1.0
        mask_inter = 1.0 - mask_intra
        
        # 내부 결합 인력(+), 그룹 간 결합 척력(-) 적용
        K_active = K_raw * mask_intra * 1.2 - K_raw * mask_inter * 0.15 + (mask_intra * 0.3 - mask_inter * 0.05)
        return K_active

    def step_euler_maruyama(self, theta, omega):
        """
        1차 Euler-Maruyama 이산 적분 1회 수행
        """
        sigma = self.get_physics_parameters()
        K_active = self._get_coupling_matrix()
        
        # 상호작용 힘 계산
        theta_diff = theta.unsqueeze(1) - theta.unsqueeze(2)  # [B, N, N]
        interaction = torch.sum(K_active * torch.sin(theta_diff), dim=-1) / self.N # [B, N]
        
        # 1차 Deterministic 위상 도함수: dtheta_dt = omega + interaction
        dtheta_dt = omega + interaction
        
        # Langevin Stochastic Step (dt의 제곱근에 비례하는 노이즈 가산)
        noise = sigma * torch.randn_like(theta) * np.sqrt(self.dt)
        theta_next = theta + dtheta_dt * self.dt + noise
        
        return theta_next

    def simulate(self, theta_0, omega, steps=48):
        """
        다단계 시간 이산 시뮬레이션 수행 및 질서 매개변수 r 출력
        """
        batch_size = theta_0.size(0)
        theta_traj = torch.zeros(batch_size, steps, self.N, device=theta_0.device)
        r_history = torch.zeros(batch_size, steps, device=theta_0.device)
        
        theta = theta_0.clone()
        
        for step in range(steps):
            theta = self.step_euler_maruyama(theta, omega)
            theta_traj[:, step] = theta
            
            # 오더 파라미터 r 계산
            r = torch.abs(torch.mean(torch.exp(1j * (theta % (2 * np.pi))), dim=-1))
            r_history[:, step] = r
            
        return theta_traj, r_history

# 사용 데모 스크립트
if __name__ == "__main__":
    batch_size = 2
    solver = StochasticKuramoto1stSolver(num_oscillators=64, dt=0.02, sigma=0.06)
    
    # 초기 상태 생성
    theta_init = torch.rand(batch_size, 64) * 2 * np.pi
    omega_init = torch.randn(batch_size, 64) * 0.15
    
    # 48일(스텝) 시뮬레이션
    traj_theta, traj_r = solver.simulate(theta_init, omega_init, steps=48)
    
    print("1차 SDE 오실레이터 시뮬레이션 완료:")
    print("위상 궤적 텐서 모양 (Theta Trajectory shape):", traj_theta.shape)  # [2, 48, 64]
    print("오더 파라미터 변동 추이 (최종 5 스텝):", traj_r[0, -5:].detach().numpy())
```
