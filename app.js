/**
 * Failure Theory Calculator — app.js
 * 3D Stress Analysis · Multi-Theory Comparison
 * Supports: von Mises, Tresca, Ductile Coulomb-Mohr,
 *           Maximum Normal Stress, Brittle Coulomb-Mohr
 */

(function () {
    'use strict';

    // ================================================================
    //  DOM REFERENCES
    // ================================================================

    const btnCalculate = document.getElementById('btn-calculate');
    const btnReset = document.getElementById('btn-reset');
    const resultsSection = document.getElementById('results-section');

    // Material inputs
    const matInputs = {
        syt: document.getElementById('input-syt'),
        syc: document.getElementById('input-syc'),
        sut: document.getElementById('input-sut'),
        suc: document.getElementById('input-suc'),
        targetSF: document.getElementById('input-target-sf'),
    };

    // Stress inputs
    const stressInputs = {
        sx: document.getElementById('input-sx'),
        sy: document.getElementById('input-sy-stress'),
        sz: document.getElementById('input-sz'),
        txy: document.getElementById('input-txy'),
        tyz: document.getElementById('input-tyz'),
        txz: document.getElementById('input-txz'),
    };

    // Error spans (only for material inputs that require validation)
    const matErrors = {
        syt: document.getElementById('error-syt'),
        syc: document.getElementById('error-syc'),
        sut: document.getElementById('error-sut'),
        suc: document.getElementById('error-suc'),
        targetSF: document.getElementById('error-target-sf'),
    };

    // Result elements
    const tensorCells = {
        xx: document.getElementById('t-xx'),
        xy: document.getElementById('t-xy'),
        xz: document.getElementById('t-xz'),
        yx: document.getElementById('t-yx'),
        yy: document.getElementById('t-yy'),
        yz: document.getElementById('t-yz'),
        zx: document.getElementById('t-zx'),
        zy: document.getElementById('t-zy'),
        zz: document.getElementById('t-zz'),
    };

    const principalEls = {
        s1: document.getElementById('res-s1'),
        s2: document.getElementById('res-s2'),
        s3: document.getElementById('res-s3'),
        tmax: document.getElementById('res-tmax'),
    };

    const theoryTbody = document.getElementById('theory-tbody');
    const envelopeCanvas = document.getElementById('envelope-chart');

    let chartInstance = null;

    // ================================================================
    //  VALIDATION
    // ================================================================

    function validateMaterial() {
        let valid = true;

        // Clear
        Object.values(matErrors).forEach(e => (e.textContent = ''));
        Object.values(matInputs).forEach(i => i.classList.remove('input-invalid'));

        // Required + positive checks
        const checks = [
            { key: 'syt', name: 'Tensile yield strength', mustBePositive: true },
            { key: 'syc', name: 'Compressive yield strength', mustBePositive: true },
            { key: 'sut', name: 'Ultimate tensile strength', mustBePositive: true },
            { key: 'suc', name: 'Ultimate compressive strength', mustBePositive: true },
            { key: 'targetSF', name: 'Target safety factor', mustBePositive: true },
        ];

        for (const chk of checks) {
            const val = matInputs[chk.key].value.trim();
            if (val === '' || isNaN(parseFloat(val))) {
                matErrors[chk.key].textContent = 'Required.';
                matInputs[chk.key].classList.add('input-invalid');
                valid = false;
            } else if (chk.mustBePositive && parseFloat(val) <= 0) {
                matErrors[chk.key].textContent = 'Must be > 0.';
                matInputs[chk.key].classList.add('input-invalid');
                valid = false;
            }
        }

        return valid;
    }

    // ================================================================
    //  PRINCIPAL STRESS SOLVER (3x3 Symmetric Eigenvalues)
    // ================================================================

    /**
     * Analytically solve the characteristic equation of a 3x3 symmetric
     * stress tensor for its three real eigenvalues (principal stresses).
     *
     * Characteristic equation: λ³ - I₁λ² + I₂λ - I₃ = 0
     *
     * Uses the trigonometric method for three guaranteed real roots.
     */
    function computePrincipalStresses(sx, sy, sz, txy, tyz, txz) {
        // Stress invariants
        const I1 = sx + sy + sz;
        const I2 = sx * sy + sy * sz + sx * sz - txy * txy - tyz * tyz - txz * txz;
        const I3 =
            sx * sy * sz +
            2 * txy * tyz * txz -
            sx * tyz * tyz -
            sy * txz * txz -
            sz * txy * txy;

        // Reduced cubic: t³ + pt + q = 0  via  λ = t + I₁/3
        const p = I1 / 3;
        const s2 = (I1 * I1 - 3 * I2) / 9;

        // Handle near-zero case (hydrostatic stress)
        if (s2 < 1e-20) {
            return [p, p, p];
        }

        const q2 = (2 * I1 * I1 * I1 - 9 * I1 * I2 + 27 * I3) / 54;
        const r = Math.sqrt(s2 * s2 * s2);
        const cosArg = Math.max(-1, Math.min(1, q2 / r));
        const phi = Math.acos(cosArg);
        const m = 2 * Math.sqrt(s2);

        const s1_val = p + m * Math.cos(phi / 3);
        const s2_val = p + m * Math.cos((phi - 2 * Math.PI) / 3);
        const s3_val = p + m * Math.cos((phi + 2 * Math.PI) / 3);

        // Sort descending: σ₁ ≥ σ₂ ≥ σ₃
        return [s1_val, s2_val, s3_val].sort((a, b) => b - a);
    }

    // ================================================================
    //  FAILURE THEORIES
    // ================================================================

    /**
     * Von Mises (Distortion Energy) — Ductile
     * σ' = √(½[(σ₁−σ₂)² + (σ₂−σ₃)² + (σ₃−σ₁)²])
     */
    function vonMises(s1, s2, s3, Syt) {
        const sigmaPrime = Math.sqrt(
            0.5 * (
                Math.pow(s1 - s2, 2) +
                Math.pow(s2 - s3, 2) +
                Math.pow(s3 - s1, 2)
            )
        );
        const fos = sigmaPrime < 1e-15 ? Infinity : Syt / sigmaPrime;
        return { effectiveStress: sigmaPrime, fos };
    }

    /**
     * Tresca (Maximum Shear Stress) — Ductile
     * τ_max = (σ_max − σ_min) / 2  →  n = Sy / (σ_max − σ_min)
     */
    function tresca(s1, s2, s3, Syt) {
        const diff = s1 - s3; // s1 is max, s3 is min after sorting
        const fos = diff < 1e-15 ? Infinity : Syt / diff;
        return { effectiveStress: diff, fos };
    }

    /**
     * Ductile Coulomb-Mohr — Ductile
     * Uses Syt (tensile yield) and Syc (compressive yield)
     */
    function ductileCoulombMohr(s1, s3, Syt, Syc) {
        let fos, effectiveStress;

        if (s1 >= 0 && s3 >= 0) {
            // All tensile
            effectiveStress = s1;
            fos = s1 < 1e-15 ? Infinity : Syt / s1;
        } else if (s1 <= 0 && s3 <= 0) {
            // All compressive
            effectiveStress = Math.abs(s3);
            fos = Math.abs(s3) < 1e-15 ? Infinity : Syc / Math.abs(s3);
        } else {
            // Mixed: σ₁ > 0, σ₃ < 0
            // 1/n = σ₁/Syt − σ₃/Syc  (note σ₃ is negative)
            const invN = s1 / Syt - s3 / Syc;
            fos = invN < 1e-15 ? Infinity : 1 / invN;
            effectiveStress = 1 / invN * Syt; // Equivalent effective stress
        }

        return { effectiveStress, fos };
    }

    /**
     * Maximum Normal Stress — Brittle
     * n = min(Sut/σ₁ if σ₁>0, Suc/|σ₃| if σ₃<0)
     */
    function maxNormalStress(s1, s3, Sut, Suc) {
        const candidates = [];
        let effectiveStress = 0;

        if (s1 > 1e-15) {
            candidates.push(Sut / s1);
            effectiveStress = s1;
        }
        if (s3 < -1e-15) {
            candidates.push(Suc / Math.abs(s3));
            if (Math.abs(s3) > effectiveStress) effectiveStress = Math.abs(s3);
        }

        const fos = candidates.length === 0 ? Infinity : Math.min(...candidates);
        return { effectiveStress, fos };
    }

    /**
     * Brittle Coulomb-Mohr — Brittle
     * Same structure as Ductile Coulomb-Mohr but uses Sut, Suc
     */
    function brittleCoulombMohr(s1, s3, Sut, Suc) {
        let fos, effectiveStress;

        if (s1 >= 0 && s3 >= 0) {
            effectiveStress = s1;
            fos = s1 < 1e-15 ? Infinity : Sut / s1;
        } else if (s1 <= 0 && s3 <= 0) {
            effectiveStress = Math.abs(s3);
            fos = Math.abs(s3) < 1e-15 ? Infinity : Suc / Math.abs(s3);
        } else {
            const invN = s1 / Sut - s3 / Suc;
            fos = invN < 1e-15 ? Infinity : 1 / invN;
            effectiveStress = 1 / invN * Sut;
        }

        return { effectiveStress, fos };
    }

    // ================================================================
    //  DISPLAY HELPERS
    // ================================================================

    function fmt(val, decimals) {
        decimals = decimals !== undefined ? decimals : 4;
        if (!isFinite(val)) return '∞';
        return val.toFixed(decimals);
    }

    function displayTensor(sx, sy, sz, txy, tyz, txz) {
        tensorCells.xx.textContent = fmt(sx, 2);
        tensorCells.xy.textContent = fmt(txy, 2);
        tensorCells.xz.textContent = fmt(txz, 2);
        tensorCells.yx.textContent = fmt(txy, 2);
        tensorCells.yy.textContent = fmt(sy, 2);
        tensorCells.yz.textContent = fmt(tyz, 2);
        tensorCells.zx.textContent = fmt(txz, 2);
        tensorCells.zy.textContent = fmt(tyz, 2);
        tensorCells.zz.textContent = fmt(sz, 2);

        // Highlight diagonal and shear
        tensorCells.xx.className = 'highlight-diag';
        tensorCells.yy.className = 'highlight-diag';
        tensorCells.zz.className = 'highlight-diag';
        tensorCells.xy.className = 'highlight-shear';
        tensorCells.xz.className = 'highlight-shear';
        tensorCells.yx.className = 'highlight-shear';
        tensorCells.yz.className = 'highlight-shear';
        tensorCells.zx.className = 'highlight-shear';
        tensorCells.zy.className = 'highlight-shear';
    }

    function displayPrincipal(s1, s2, s3) {
        principalEls.s1.textContent = fmt(s1, 4);
        principalEls.s2.textContent = fmt(s2, 4);
        principalEls.s3.textContent = fmt(s3, 4);
        principalEls.tmax.textContent = fmt((s1 - s3) / 2, 4);
    }

    function displayTheoryTable(theories, targetSF) {
        theoryTbody.innerHTML = '';

        theories.forEach(t => {
            const row = document.createElement('tr');
            const isSafe = t.fos >= 1;
            const meetsTarget = t.fos >= targetSF;

            row.className = isSafe ? 'row-safe' : 'row-fail';

            // Interpretation
            let interpText, interpClass;
            if (!isFinite(t.fos)) {
                interpText = 'No stress applied';
                interpClass = 'meets';
            } else if (meetsTarget) {
                interpText = `✓ Meets target (n ≥ ${fmt(targetSF, 2)})`;
                interpClass = 'meets';
            } else if (isSafe) {
                interpText = `⚠ Below target (n < ${fmt(targetSF, 2)})`;
                interpClass = 'marginal';
            } else {
                interpText = '✗ Predicts failure';
                interpClass = 'below';
            }

            row.innerHTML = `
                <td><span class="theory-name">${t.name}</span></td>
                <td><span class="behavior-badge ${t.behavior}">${t.behavior}</span></td>
                <td><span class="effective-stress-value">${fmt(t.effectiveStress, 4)}</span></td>
                <td><span class="fos-value" style="color: ${isSafe ? 'var(--safe-green)' : 'var(--danger-red)'}">${fmt(t.fos, 4)}</span></td>
                <td><span class="status-pill ${isSafe ? 'safe' : 'failure'}">${isSafe ? '✓ Safe' : '✗ Failure'}</span></td>
                <td><span class="interpretation-text ${interpClass}">${interpText}</span></td>
            `;

            theoryTbody.appendChild(row);
        });
    }

    // ================================================================
    //  FAILURE ENVELOPE CHART (Chart.js)
    // ================================================================

    function generateVonMisesEllipse(Sy, numPoints) {
        numPoints = numPoints || 200;
        const points = [];
        for (let i = 0; i <= numPoints; i++) {
            const t = (2 * Math.PI * i) / numPoints;
            const s1 = Sy * (Math.cos(t) + Math.sin(t) / Math.sqrt(3));
            const s2 = Sy * (Math.cos(t) - Math.sin(t) / Math.sqrt(3));
            points.push({ x: s1, y: s2 });
        }
        return points;
    }

    function generateTrescaHexagon(Sy) {
        return [
            { x: Sy, y: 0 },
            { x: Sy, y: Sy },
            { x: 0, y: Sy },
            { x: -Sy, y: 0 },
            { x: -Sy, y: -Sy },
            { x: 0, y: -Sy },
            { x: Sy, y: 0 },    // close the hexagon
        ];
    }

    function drawEnvelopeChart(Syt, sigma1, sigma2) {
        const vmPoints = generateVonMisesEllipse(Syt);
        const trescaPoints = generateTrescaHexagon(Syt);

        // Determine axis range
        const maxStress = Math.max(
            Syt * 1.35,
            Math.abs(sigma1) * 1.3,
            Math.abs(sigma2) * 1.3,
            100
        );
        const axisMin = -maxStress;
        const axisMax = maxStress;

        // Axis reference lines
        const xAxisLine = [{ x: axisMin, y: 0 }, { x: axisMax, y: 0 }];
        const yAxisLine = [{ x: 0, y: axisMin }, { x: 0, y: axisMax }];

        const datasets = [
            // X axis reference
            {
                label: 'X axis',
                data: xAxisLine,
                showLine: true,
                fill: false,
                borderColor: 'rgba(0,0,0,0.08)',
                borderWidth: 1,
                pointRadius: 0,
                borderDash: [],
                order: 5,
            },
            // Y axis reference
            {
                label: 'Y axis',
                data: yAxisLine,
                showLine: true,
                fill: false,
                borderColor: 'rgba(0,0,0,0.08)',
                borderWidth: 1,
                pointRadius: 0,
                borderDash: [],
                order: 5,
            },
            // Von Mises envelope
            {
                label: 'von Mises Envelope',
                data: vmPoints,
                showLine: true,
                fill: false,
                borderColor: '#003366',
                backgroundColor: 'rgba(0, 51, 102, 0.05)',
                borderWidth: 2.5,
                pointRadius: 0,
                tension: 0,
                order: 2,
                fill: true,
            },
            // Tresca hexagon
            {
                label: 'Tresca Hexagon',
                data: trescaPoints,
                showLine: true,
                fill: false,
                borderColor: '#c8a951',
                borderWidth: 2.5,
                pointRadius: 0,
                borderDash: [10, 5],
                tension: 0,
                order: 3,
            },
            // Operating point
            {
                label: `Operating Point (${fmt(sigma1, 1)}, ${fmt(sigma2, 1)})`,
                data: [{ x: sigma1, y: sigma2 }],
                pointRadius: 10,
                pointHoverRadius: 13,
                pointBackgroundColor: '#e74c3c',
                pointBorderColor: '#ffffff',
                pointBorderWidth: 3,
                showLine: false,
                order: 1,
            },
        ];

        if (chartInstance) {
            chartInstance.destroy();
        }

        chartInstance = new Chart(envelopeCanvas, {
            type: 'scatter',
            data: { datasets },
            options: {
                responsive: true,
                maintainAspectRatio: true,
                aspectRatio: 1,
                animation: {
                    duration: 800,
                    easing: 'easeOutQuart',
                },
                plugins: {
                    legend: {
                        position: 'top',
                        labels: {
                            usePointStyle: true,
                            padding: 16,
                            font: { family: "'Inter', sans-serif", size: 12, weight: '500' },
                            filter: function (item) {
                                // Hide axis reference lines from legend
                                return item.text !== 'X axis' && item.text !== 'Y axis';
                            },
                        },
                    },
                    title: {
                        display: true,
                        text: 'Failure Envelope — Plane Stress (σ₃ = 0)',
                        font: { family: "'Inter', sans-serif", size: 15, weight: '700' },
                        color: '#003366',
                        padding: { bottom: 16 },
                    },
                    tooltip: {
                        backgroundColor: '#00264d',
                        titleFont: { family: "'Inter', sans-serif", weight: '600' },
                        bodyFont: { family: "'JetBrains Mono', monospace", size: 12 },
                        padding: 10,
                        cornerRadius: 8,
                        callbacks: {
                            label: function (ctx) {
                                return `σ₁ = ${ctx.parsed.x.toFixed(2)}, σ₂ = ${ctx.parsed.y.toFixed(2)} MPa`;
                            },
                        },
                    },
                },
                scales: {
                    x: {
                        type: 'linear',
                        min: axisMin,
                        max: axisMax,
                        title: {
                            display: true,
                            text: 'σ₁ (MPa)',
                            font: { family: "'Inter', sans-serif", size: 13, weight: '600' },
                            color: '#003366',
                        },
                        grid: { color: 'rgba(0,0,0,0.05)' },
                        ticks: {
                            font: { family: "'JetBrains Mono', monospace", size: 10 },
                            color: '#7a8494',
                        },
                    },
                    y: {
                        type: 'linear',
                        min: axisMin,
                        max: axisMax,
                        title: {
                            display: true,
                            text: 'σ₂ (MPa)',
                            font: { family: "'Inter', sans-serif", size: 13, weight: '600' },
                            color: '#003366',
                        },
                        grid: { color: 'rgba(0,0,0,0.05)' },
                        ticks: {
                            font: { family: "'JetBrains Mono', monospace", size: 10 },
                            color: '#7a8494',
                        },
                    },
                },
            },
        });
    }

    // ================================================================
    //  MAIN CALCULATION HANDLER
    // ================================================================

    function runCalculation() {
        if (!validateMaterial()) return;

        // Read material properties
        const Syt = parseFloat(matInputs.syt.value);
        const Syc = parseFloat(matInputs.syc.value);
        const Sut = parseFloat(matInputs.sut.value);
        const Suc = parseFloat(matInputs.suc.value);
        const targetSF = parseFloat(matInputs.targetSF.value);

        // Read stresses (default to 0)
        const sx = parseFloat(stressInputs.sx.value) || 0;
        const sy = parseFloat(stressInputs.sy.value) || 0;
        const sz = parseFloat(stressInputs.sz.value) || 0;
        const txy = parseFloat(stressInputs.txy.value) || 0;
        const tyz = parseFloat(stressInputs.tyz.value) || 0;
        const txz = parseFloat(stressInputs.txz.value) || 0;

        // --- Display tensor ---
        displayTensor(sx, sy, sz, txy, tyz, txz);

        // --- Compute principal stresses ---
        const [s1, s2, s3] = computePrincipalStresses(sx, sy, sz, txy, tyz, txz);
        displayPrincipal(s1, s2, s3);

        // --- Compute all theories ---
        const vmResult = vonMises(s1, s2, s3, Syt);
        const trResult = tresca(s1, s2, s3, Syt);
        const dcmResult = ductileCoulombMohr(s1, s3, Syt, Syc);
        const mnsResult = maxNormalStress(s1, s3, Sut, Suc);
        const bcmResult = brittleCoulombMohr(s1, s3, Sut, Suc);

        const theories = [
            {
                name: 'Distortion Energy (von Mises)',
                behavior: 'ductile',
                effectiveStress: vmResult.effectiveStress,
                fos: vmResult.fos,
            },
            {
                name: 'Max Shear Stress (Tresca)',
                behavior: 'ductile',
                effectiveStress: trResult.effectiveStress,
                fos: trResult.fos,
            },
            {
                name: 'Ductile Coulomb-Mohr',
                behavior: 'ductile',
                effectiveStress: dcmResult.effectiveStress,
                fos: dcmResult.fos,
            },
            {
                name: 'Maximum Normal Stress',
                behavior: 'brittle',
                effectiveStress: mnsResult.effectiveStress,
                fos: mnsResult.fos,
            },
            {
                name: 'Brittle Coulomb-Mohr',
                behavior: 'brittle',
                effectiveStress: bcmResult.effectiveStress,
                fos: bcmResult.fos,
            },
        ];

        displayTheoryTable(theories, targetSF);

        // --- Draw failure envelope ---
        drawEnvelopeChart(Syt, s1, s2);

        // --- Show results ---
        resultsSection.classList.remove('hidden');
        resultsSection.style.animation = 'none';
        void resultsSection.offsetHeight;
        resultsSection.style.animation = '';

        // Scroll to results
        setTimeout(function () {
            resultsSection.scrollIntoView({ behavior: 'smooth', block: 'start' });
        }, 100);
    }

    // ================================================================
    //  RESET
    // ================================================================

    function resetAll() {
        // Clear material inputs
        Object.values(matInputs).forEach(i => {
            if (i.id === 'input-target-sf') {
                i.value = '1';
            } else {
                i.value = '';
            }
            i.classList.remove('input-invalid');
        });
        Object.values(matErrors).forEach(e => (e.textContent = ''));

        // Reset stress inputs to 0
        Object.values(stressInputs).forEach(i => (i.value = '0'));

        // Hide results
        resultsSection.classList.add('hidden');

        // Destroy chart
        if (chartInstance) {
            chartInstance.destroy();
            chartInstance = null;
        }
    }

    // ================================================================
    //  EVENT LISTENERS
    // ================================================================

    btnCalculate.addEventListener('click', runCalculation);
    btnReset.addEventListener('click', resetAll);

    // Allow Enter key to trigger calculation from any input
    document.addEventListener('keydown', function (e) {
        if (e.key === 'Enter' && e.target.tagName === 'INPUT') {
            e.preventDefault();
            runCalculation();
        }
    });

})();
