#[cfg(test)]
mod tests;

// Port of:
// McCoy, R. L. (1981).
// " MC DRAG"-a computer program for estimating the drag coefficients of projectiles.
// ARMY BALLISTIC RESEARCH LAB ABERDEEN PROVING GROUND MD.

type Float = f64;

#[derive(Debug, Clone, Copy)]
pub enum BoundaryLayerCode
{
    LL,
    LT,
    TT,
}

#[derive(Debug, Clone)]
pub struct MCDRAG
// where Float
{
    pub d1 : Float, // Projectile reference diameter (mm)
    pub l1 : Float, // Projectile length (calibers)
    pub l2 : Float, // Nose length (calibers)
    pub r1 : Float, // RT/R (headshape parameter)
    pub l3 : Float, // Boattail length (calibers)
    pub d2 : Float, // Base diameter (calibers)
    pub d3 : Float, // Meplat diameter (calibers)
    pub d4 : Float, // Rotating band diameter (calibers)
    pub x1 : Float, // Center of gravity location (optional, calibers from nose)
    pub k :  BoundaryLayerCode,
}

impl MCDRAG
{
    pub fn new(
        projectile_reference_diameter_mm : Float,
        projectile_length_calibers : Float,
        nose_length_calibers : Float,
        rt_r_headshape_parameter : Float,
        boattail_length_calibers : Float,
        base_diameter_calibers : Float,
        meplat_diameter_calibers : Float,
        rotating_band_diameter : Float,
        center_of_gravity_location_calibers_from_nose : Float,
        boundary_layer_code : BoundaryLayerCode,
    ) -> Self
    {
        Self {
            d1 : projectile_reference_diameter_mm,
            l1 : projectile_length_calibers,
            l2 : nose_length_calibers,
            r1 : rt_r_headshape_parameter,
            l3 : boattail_length_calibers,
            d2 : base_diameter_calibers,
            d3 : meplat_diameter_calibers,
            d4 : rotating_band_diameter,
            x1 : center_of_gravity_location_calibers_from_nose,
            k :  boundary_layer_code,
        }
    }

    // Generates the drag coefficient given a velocity (in machs) through any fluid
    pub fn gen_unclamped(&self, mach : Float) -> Float
    {
        let c17 : Float;
        let c11 : Float;
        let c12 : Float;
        let c13 : Float;
        let c14 : Float;
        let r4 : Float;
        let t1 : Float;
        let m2 : Float;
        let r2 : Float;
        let r3 : Float;
        let c7 : Float;
        let c8 : Float;
        let d5 : Float;
        let s1 : Float;
        let s2 : Float;
        let s3 : Float;
        let c9 : Float;
        let c10 : Float;
        let b2 : Float;
        let b : Float;
        let mut z : Float;
        let s4 : Float;
        let c15 : Float;
        let p5 : Float;
        let c16 : Float;
        let c18 : Float;
        let t2 : Float;
        let t3 : Float;
        let e1 : Float;
        let b4 : Float;
        let b3 : Float;
        let a12 : Float;
        let a11 : Float;
        let e2 : Float;
        let x3 : Float;
        let a1 : Float;
        let r5 : Float;
        let e3 : Float;
        let a2 : Float;
        let p2 : Float;
        let p4 : Float;
        let x2 : Float;

        // These were initially arrays storing values at points at velocity (like mach)
        let c1 : Float;
        let c2 : Float;
        let c3 : Float;
        let c4 : Float;
        let c5 : Float;
        let c6 : Float;
        let p1 : Float;

        // Line 74
        t1 = (1.0 - self.d3) / self.l2;
        m2 = mach * mach;
        r2 = 23296.3 * mach * self.l1 * self.d1;
        r3 = 0.4343 * r2.ln();
        c7 = (1.328 / r2.sqrt()) * ((1.0 + 0.12 * m2).powf(-0.12));
        c8 = (0.455 / r3.powf(2.58)) * ((1.0 + 0.21 * m2).powf(-0.32));
        d5 = 1.0 + ((0.333 + (0.02 / (self.l2 * self.l2))) * self.r1);
        s1 = 1.5708 * self.l2 * d5 * (1.0 + (1.0 / (8.0 * self.l2 * self.l2)));
        s2 = 3.1416 * (self.l1 - self.l2);
        s3 = s1 + s2;

        // Line 84
        match self.k
        {
            // Line 87
            BoundaryLayerCode::LL =>
            {
                c9 = 1.2732 * s3 * c7;
                c10 = c9;
            },
            // Line 90
            BoundaryLayerCode::LT =>
            {
                c9 = 1.2732 * s3 * c7;
                c10 = 1.2732 * s3 * c8;
            },
            // Line 93
            BoundaryLayerCode::TT =>
            {
                c9 = 1.2732 * s3 * c8;
                c10 = c9;
            },
        };

        c3 = (c9 * s1 + c10 * s2) / s3;
        c15 = (m2 - 1.0) / (2.4 * m2);

        if mach > 1.0
        {
            p5 = (1.2 * m2).powf(3.5) * (6.0 / (7.0 * m2 - 1.0)).powf(2.5);
        }
        else
        {
            p5 = (1.0 + 0.2 * m2).powf(3.5);
        }

        // Line 102
        c16 = (1.122 * (p5 - 1.0) * self.d3).powi(2) / m2;

        if mach <= 0.91
        {
            c18 = 0.0;
        }
        else if mach >= 1.41
        {
            c18 = 0.85 * c16;
        }
        else
        {
            c18 = (0.254 + 2.88 * c15) * c16;
        }

        // Line 110
        if mach < 1.0
        {
            p2 = 1.0 / (1.0 + 0.1875 * m2 + 0.0531 * (m2).powi(2));
        }
        else
        {
            p2 = 1.0 / (1.0 + 0.2477 * m2 + 0.0345 * (m2).powi(2));
        }

        p4 = (1.0 + 9.000001e-02 * m2 * (1.0 - (self.l2 - self.l1).exp())) * (1.0 + 0.25 * m2 * (1.0 - self.d2));

        // Line 116
        p1 = p2 * p4;
        if p1 >= 0.0
        {
            c6 = (1.4286 * (1.0 - p1) * (self.d2 * self.d2)) / m2;
        }
        else
        {
            c6 = 0.0;
        }

        if mach < 0.95
        {
            c4 = mach.powf(12.5) * (self.d4 - 1.0);
        }
        else
        {
            c4 = (0.21 + 0.28 / m2) * (self.d4 - 1.0);
        }

        // Line 126
        if mach > 1.0
        {
            // line 146
            // supersonic speeds
            b2 = m2 - 1.0;
            b = b2.sqrt();
            z = b;
            s4 = 1.0 + 0.368 * t1.powf(1.85);

            // line 150
            if mach < s4
            {
                z = (s4 * s4 - 1.0).sqrt();
            }

            c11 = 0.7156 - 0.5313 * self.r1 + 0.595 * self.r1 * self.r1;
            c12 = 0.0796 + 0.0779 * self.r1;
            c13 = 1.587 + 0.049 * self.r1;
            c14 = 0.1122 + 0.1658 * self.r1;
            r4 = 1.0 / (z * z);
            c17 = (c11 - c12 * (t1 * t1)) * r4 * (t1 * z).powf(c13 + c14 * t1);
            c2 = c17 + c18;

            // line 159
            if self.l3 <= 0.0
            {
                c5 = 0.0;
            }
            else
            {
                t2 = (1.0 - self.d2) / (2.0 * self.l3);

                if mach <= 1.1
                {
                    t3 = 2.0 * t2 * t2 + t2 * t2 * t2;
                    e1 = (-2.0 * self.l3).exp();
                    b4 = 1.0 - e1 + 2.0 * t2 * ((e1 * (self.l3 + 0.5)) - 0.5);
                    c5 = 2.0 * t3 * b4 * (1.774 - 9.3 * c15);
                }
                else
                {
                    b3 = 0.85 / b;
                    a12 = (5.0 * t1) / (6.0 * b) + (0.5 * t1) * (0.5 * t1) - (0.7435 / m2) * (t1 * mach).powf(1.6);
                    a11 = (1.0 - ((0.6 * self.r1) / mach)) * a12;
                    e2 = (((-1.1952) / mach) * (self.l1 - self.l2 - self.l3)).exp();
                    x3 = ((2.4 * m2 * m2 - 4.0 * b2) * (t2 * t2)) / (2.0 * b2 * b2);
                    a1 = a11 * e2 - x3 + ((2.0 * t2) / b);
                    r5 = 1.0 / b3;
                    e3 = (-b3 * self.l3).exp();
                    a2 = 1.0 - e3 + (2.0 * t2 * (e3 * (self.l3 + r5) - r5));
                    c5 = 4.0 * a1 * t2 * a2 * r5;
                }
            }
        }
        else
        {
            if self.l3 <= 0.0 || mach <= 0.85
            {
                // line 129
                c5 = 0.0;
            }
            else
            {
                // line 131 (moved check into previous if)

                t2 = (1.0 - self.d2) / (2.0 * self.l3);
                t3 = 2.0 * t2 * t2 + t2 * t2 * t2;
                e1 = (-2.0 * self.l3).exp();
                b4 = 1.0 - e1 + 2.0 * t2 * ((e1 * (self.l3 + 0.5)) - 0.5);
                c5 = 2.0 * t3 * b4 * (1.0 / (0.564 + 1250.0 * c15 * c15));
            }

            x2 = (1.0 + 0.552 * t1.powf(0.8)).powf(-0.5);
            // line 138
            if mach <= x2
            {
                c2 = c18; // should be 'c2[i] = c17 + c18', but c17 is already '0'.

            // goto 181
            }
            else
            {
                c17 = 0.368 * t1.powf(1.8) + 1.6 * t1 * c15;
                c2 = c17 + c18;

                // goto 181
            }
        }

        // line 181
        // c1 === cdo (refer to page.10)
        c1 =
            Self::fix_float(c2) + Self::fix_float(c3) + Self::fix_float(c4) + Self::fix_float(c5) + Self::fix_float(c6);

        c1
    }

    pub fn gen(&self, mach : Float) -> Float { self.gen_unclamped(mach.clamp(0.5, 5.0)) }

    #[inline]
    fn fix_float(v : Float) -> Float
    {
        if v.is_infinite() || v.is_nan()
        {
            0.0
        }
        else
        {
            v
        }
    }
}
