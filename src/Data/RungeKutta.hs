module Data.RungeKutta where

import Linear

data ButcherTableau f g c = ButcherTableau { w :: f (g c)
                                           , c :: f c
                                           , b :: g c
                                           }
    deriving (Show)


butcherTableau11 = ButcherTableau { w = V1 (V1 1)
                                  , c = V1 0
                                  , b = V1 1
                                  }
-- lack of stability and accuracy, only serves as a simple example
eulersMethod = butcherTableau11

-- alpha \in (0, 1]
butcherTableau22 alpha = ButcherTableau { w = V2 (V2 0 0) (V2 alpha 0)
                                        , c = V2 0 alpha
                                        , b = V2 (1-alpha') alpha'
                                        }
    where alpha' = 1/(2*alpha)

-- explicit midpoint method
midpointRule = butcherTableau22 $ 1/2

-- explicit trapezoid rule
heunsMethod = butcherTableau22 1

-- has a minimum local error bound
ralstonsMethod = butcherTableau22 $ 2/3

-- u, v \in (0, 1]
butcherTableau33 u v = ButcherTableau { w = V3 (V3 0 0 0) (V3 u 0 0) (V3 (v-uv') uv' 0)
                                      , c = V3 0 u v
                                      , b = V3 (1-u'-v') u' v'
                                      }
    where uv' = v*(v-u)/(u*(2-3*u))
          u' = (2-3*v)/(6*u*(u-v))
          v' = (2-3*u)/(6*v*(v-u))

kuttas3Method = butcherTableau33 (1/2) 1

-- The "original" Runge-Kutta method
kuttas4Method = ButcherTableau { w = V4 (V4 0 0 0 0) (V4 (1/2) 0 0 0) (V4 0 (1/2) 0 0) (V4 0 0 1 0)
                               , c = V4 0 (1/2) (1/2) 1
                               , b = V4 (1/6) (1/3) (1/3) (1/6)
                               }



-- dx/dt = f(x, t)
-- x(t0) = x0
-- t \in T = [t0, tH]

-- generalize [a] to (Applicative f) => a -- or something
rk :: (Num t) => ButcherTableau m n v -> t -> (v -> t -> v) -> (t, v) -> [(t, v)]
rk bt h f (x0, t0) = (x0, t0) : rk' bt h f (x0, t0)
    where
        rk' bt@ButcherTableau{ w=w, c=c, b=b} h f (t, x) = (t', x') : rk' bt h f (t', x')
            where
                t' = t + h
                x' = x + h * dot b y
                where
                    y :: f v
                    y = f <$> (x + h * x_) <*> (t + h * c)
                    x_ :: f v
                    x_ = .. recursive since y_j -- x_i = sum^{i-1}_{j=1} w_ij y_j
