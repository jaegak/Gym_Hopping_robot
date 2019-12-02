"""
개인적인 GYM 환경 구성을 해보기 위하여 1984년도 MARC H. RAIBERT의 Hopping in Legged Systems - Modeling and Simulation for the Two-Dimenional One-Legged Case 논문에서 나오는
one-legged hopping machine을 gym 환경 상에서 구축하는 코드를 작성하고자 함
최종적인 목표는 구축한 환경을 학습을 통해 Hopping 로봇이 목표 높이로 제자리 점프하는 것과 목표 속도로 이동할 수 있는 2가지 case에 대하여 학습을 하고자 한다.
코드 작성에 참고한 자료는 gym 내에 예제로 올라와있는 cartpole.py를 참고하였다.
또한 이전에 matlab을 통해 simulation 구현시 사용한 값들을 이용하였다.

action은 input인 Torque와 cai의 변화 값으로 주었다.
"""

import math
import gym
from gym import spaces, logger
from gym.utils import seeding
import numpy as np
import numpy.linalg as lin

class HoppingRobotEnv(gym.Env):
    metadata = {
        'render.modes': ['human','rgb_array'],
        'video.frames_per_second' : 50
    }

    def __init__(self):
        self.gravity = 9.8
        self.mass_link_1 = 1
        self.mass_link_2 = 10
        self.inertia_link_1 = 1
        self.inertia_link_2 = 10
        self.length_r1 = 0.5
        self.length_r2 = 0.4

        self.k0 = 1
        self.K_L = 1000
        self.K_L2 = 100000
        self.K_G = 10000
        self.B_L2 = 125
        self.B_G = 75
        self.K_P1 = 1800
        self.K_P2 = 1200
        self.K_V1 = 200
        self.K_V2 = 60


        self.Torque_mag = 1
        self.cai_dot_mag = 0.005
        self.tau = 0.001
        self.kinematics_integrator = 'euler'

        # Angle at which to fail the episode
        self.theta1_threshold_radians = 45 * 2 * math.pi / 360
        self.theta2_threshold_radians = 12 * 2 * math.pi / 360
        self.x0_threshold = 5
        self.w_threshold_min = 0
        self.w_threshold_max = 1

        self.cai_max = 0.15
        self.cai_init = 0.075
        self.cai_min = 0

        high = np.array([
            self.x_threshold *2,
            np.finfo(np.float32).max,
            self.theta1_threshold_radians *2,
            np.finfo(np.float32).max,
            self.theta2_threshold_radians *2,
            np.finfo(np.float32).max])

        self.action_space = spaces.Discrete(2)
        self.observation_space = spaces.Box(-high, high, dtype=np.float32)

        self.seed()
        self.viewer = None
        self.state = None

        self.steps_beyond_done = None

    def seed(self, seed = None):
        self.np_random, seed = seeding.np_random(seed)
        return [seed]

    def step(self, action):
        assert self.action_space.contains(action), "%r (%s) invalid"%(action, type(action))
        state = self.state
        x0, x0_dot, y0, y0_dot, w, w_dot, theta1, theta1_dot, theta2, theta2_dot, TouchDown, cai = state
        W=w-self.length_r1

        if y0 < 0:
            TouchDown = 1
        else:
            TouchDown = 0

        F_K=self.K_L2*(self.k0-w+cai)-self.B_L2*w_dot
        if (self.k0-w+cai) > 0:
            F_K = self.K_L*(self.k0 - w + cai)

        F_x = 0
        F_y = 0
        if y0 < 0:
            F_x = -K_G*(x0-x_TouchDown)-B_G*x0_dot
            F_y = -K_G*y0-B_G*y0_dot

        # action은 0~8 로 9개로 -값, 0 , +값으로 3가지 값을 가진 2가지 변수를 구성

        if action==0|action==1|action==2:
            Torque = -self.Torque_mag
        elif action==4|action==5|action==3:
            Torque = 0
        else:
            Torque = self.Torque_mag

        if action==0|action==3|action==6:
            cai_dot = -self.cai_dot_mag
        elif action==1|action==4|action==7:
            cai_dot= 0
        else:
            cai_dot = self.cai_dot_mag



        A11 = math.cos(theta1)*(self.mass_link_2*W*w + self.inertia_link_1)
        A12 = self.mass_link_2*self.length_r2*W*math.cos(theta2)
        A13 = self.mass_link_2*W
        A14 = 0
        A15 = self.mass_link_2*W*math.sin(theta1)

        B1 = self.mass_link_2*W*(w*theta1_dot^2*math.sin(theta1) + self.length_r2*theta2_dot^2*math.sin(theta2)
                                 - w*w_dot*theta1_dot*math.cos(theta1)) - F_x*self.length_r1*(math.cos(theta1))^2\
             + F_y*self.length_r1*math.cos(theta1)*math.sin(theta1) - Torque*math.cos(theta1) + F_K*W*math.sin(theta1)

        A21 = -math.sin(theta1)*(self.mass_link_2*W*w+self.inertia_link_1)
        A22 = -self.mass_link_2*self.length_r2*W*math.sin(theta2)
        A23 = 0
        A24 = self.mass_link_2*W
        A25 = self.mass_link_2*W*math.cos(theta1)

        B2 = self.mass_link_2*W*(w*theta1_dot^2*math.cos(theta1) + self.length_r2*theta2_dot^2*math.cos(theta2)
                                 + 2*w_dot*theta1_dot*math.sin(theta1) - self.gravity) \
             + F_x*self.length_r1*math.cos(theta1)*math.sin(theta1)\
             - F_y*self.length_r1*(math.sin(theta1))^2 + Torque*math.sin(theta1) + F_K*W*math.cos(theta1)

        A31 = math.cos(theta1)*(self.mass_link_1*self.length_r1*W - self.inertia_link_1)
        A32 = 0
        A33 = self.mass_link_1*W
        A34 = 0
        A35 = 0

        B3 = W*(F_x - F_K*math.sin(theta1) + self.mass_link_1*self.length_r1*theta1_dot^2*math.sin(theta1))\
             - math.cos(theta1)*(-F_x*self.length_r1*math.cos(theta1) + F_y*self.length_r1*math.sin(theta1) - Torque)

        A41 = -math.sin(theta1)*(self.mass_link_1*self.length_r1*W-self.inertia_link_1)
        A42 = 0
        A43 = 0
        A44 = self.mass_link_1*W
        A45 = 0

        B4 = W*(F_y - F_K*math.cos(theta1) - self.mass_link_1*self.gravity
                + self.mass_link_1*self.length_r1*theta1_dot^2*math.cos(theta1))\
             + math.sin(theta1)*(-F_x*self.length_r1*math.cos(theta1) + F_y*self.length_r1*math.sin(theta1) - Torque)

        A51 = -math.cos(theta2-theta1)*self.inertia_link_1*self.length_r2
        A52 = self.inertia_link_2*W
        A53 = 0
        A54 = 0
        A55 = 0

        B5 = W*(F_K*self.length_r2*math.sin(theta2-theta1) + Torque)\
             - self.length_r2*math.cos(theta2-theta1)*(self.length_r1*F_y*math.sin(theta1) - self.length_r1*F_x*math.cos(theta1) - Torque)

        matrix_A = np.array([A11, A12, A13, A14, A15],[A21, A22, A23, A24, A25],[A31, A32, A33, A34, A35],[A41, A42, A43, A44, A45],[A51, A52, A53, A54, A55])
        matrix_B = np.array([B1],[B2],[B3],[B4],[B5])
        inverse_matrix_A = lin.inv(matrix_A)

        matrix_C = np.dot(inverse_matrix_A, matrix_B)

        theta1_acc = matrix_C[0]
        theta2_acc = matrix_C[1]
        x0_acc = matrix_C[2]
        y0_acc = matrix_C[3]
        w_acc = matrix_C[4]


        if self.kinematics_integrator == 'euler':
            x0  = x0 + self.tau * x0_dot
            x0_dot = x0_dot + self.tau * x0_acc
            y0  = y0 + self.tau * y0_dot
            y0_dot = y0_dot + self.tau * y0_acc
            w  = w + self.tau * w_dot
            w_dot = w_dot + self.tau * w_acc
            theta1 = theta1 + self.tau * theta1_dot
            theta1_dot = theta1_dot + self.tau * theta1_acc
            theta2 = theta2 + self.tau * theta2_dot
            theta2_dot = theta2_dot + self.tau * theta2_acc
            cai = cai + self.tau * cai_dot
        else: # semi-implicit euler
            x0_dot = x0_dot + self.tau * x0_acc
            x0  = x0 + self.tau * x0_dot
            y0_dot = y0_dot + self.tau * y0_acc
            y0  = y0 + self.tau * y0_dot
            w_dot = w_dot + self.tau * w_acc
            w  = w + self.tau * w_dot
            theta1_dot = theta1_dot + self.tau * theta1_acc
            theta1 = theta1 + self.tau * theta1_dot
            theta2_dot = theta2_dot + self.tau * theta2_acc
            theta2 = theta2 + self.tau * theta2_dot
            cai = cai + self.tau * cai_dot
        self.state = (x0,x0_dot,y0,y0_dot,w,w_dot,theta1,theta1_dot,theta2,theta2_dot,TouchDown,cai)

        done =  x0 < -self.x0_threshold \
                or x0 > self.x0_threshold \
                or w < self.w_threshold_min \
                or w > self.w_threshold_max \
                or theta1 < -self.theta1_threshold_radians \
                or theta1 > self.theta1_threshold_radians \
                or theta2 < -self.theta2_threshold_radians \
                or theta2 > self.theta2_threshold_radians \
                or cai < self.cai_min \
                or cai > self.cai_max

        done = bool(done)

        if not done:
            reward = 1.0
        elif self.steps_beyond_done is None:
            # Pole just fell!
            self.steps_beyond_done = 0
            reward = 1.0
        else:
            if self.steps_beyond_done == 0:
                logger.warn("You are calling 'step()' even though this environment has already returned done = True. You should always call 'reset()' once you receive 'done = True' -- any further steps are undefined behavior.")
            self.steps_beyond_done += 120
            reward = 0.0

        return np.array(self.state), reward, done, {}

    def reset(self):
        self.state = self.np_random.uniform(low=-0.05, high=0.05, size=(12,))
        self.state[10] = 0
        self.state[11] = self.cai_init

        self.steps_beyond_done = None
        return np.array(self.state)

    def close(self):
        if self.viewer:
            self.viewer.close()
            self.viewer = None



