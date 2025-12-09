/** JavaScript で使われる関数などの説明 */
// PI...円周率
// sqrt()...平方根
// pow(底, 指数)...累乗
// randomGaussian(mean, sd)...正規分布に従う乱数を生成

/** グローバル変数（シミュレーション用） */
let params = {
    omega_1: 10.0,
    t_pulse1: 4.5,
    t_pulse2: 9.0,
    T2_star: 5.0,
    T2: 100.0
};
// 初期値を保存しておく（リセット用）
const initialParams = { ...params };

const tau = 30;
let sliders = {};
let mySim;

/** ページ読み込み時の処理 */
function setup() {
    new Canvas().create();
    createInterface();
    updateSimulation();
}

/** 画面の更新処理 */
function updateSimulation() {
    background("#000000");

    // スライダー値の読み取り
    params.omega_1  = sliders.omega_1.value();
    params.t_pulse1 = sliders.t_pulse1.value();
    params.t_pulse2 = sliders.t_pulse2.value();
    params.T2_star  = sliders.T2_star.value();
    params.T2       = sliders.T2.value();

    drawLabels();

    // シミュレーター計算
    mySim = new EnsembleSimulator(
        params.omega_1, 
        params.t_pulse1, 
        params.t_pulse2, 
        params.T2_star, 
        params.T2
    );
    
    // グラフ描画
    let graph = new Graph();
    graph.create();

    const pulseMag = mySim.calculation_PulseMag();
    graph.drawPulse(pulseMag, params.t_pulse1, 0);
    graph.drawPulse(pulseMag, params.t_pulse2, tau);
    
    graph.drawFID(mySim, params.t_pulse1, tau);
    graph.drawSE(mySim, params.t_pulse1, tau);
}

function draw() { noLoop(); }

// ============================================
// UI作成・管理
// ============================================
function createInterface() {
    let adj = new Adjuster();
    // スライダー作成
    sliders.omega_1  = adj.makeSlider("omega_1", 7.5, 12.5, params.omega_1, 0.1, 40);
    sliders.t_pulse1 = adj.makeSlider("t_pulse1", 0.0, 10.0, params.t_pulse1, 0.1, 90);
    sliders.t_pulse2 = adj.makeSlider("t_pulse2", 0.0, 20.0, params.t_pulse2, 0.1, 140);
    sliders.T2_star  = adj.makeSlider("T2* (Inhomo)", 1.0, 20.0, params.T2_star, 1.0, 190);
    sliders.T2       = adj.makeSlider("T2 (True)", 10.0, 500.0, params.T2, 10.0, 240);

    // ★追加：リセットボタン
    let resetBtn = createButton('Reset');
    resetBtn.position(adj.canvasWidth - 190, 300); // 一番下に配置
    resetBtn.style('width', '160px');
    resetBtn.style('height', '30px');
    resetBtn.style('background-color', '#444');
    resetBtn.style('color', 'white');
    resetBtn.style('border', '1px solid #888');
    resetBtn.style('cursor', 'pointer');
    
    // ボタンを押した時の動作
    resetBtn.mousePressed(() => {
        sliders.omega_1.value(initialParams.omega_1);
        sliders.t_pulse1.value(initialParams.t_pulse1);
        sliders.t_pulse2.value(initialParams.t_pulse2);
        sliders.T2_star.value(initialParams.T2_star);
        sliders.T2.value(initialParams.T2);
        updateSimulation();
    });
}

function drawLabels() {
    let adj = new Adjuster();
    adj.drawLabel("omega_1", params.omega_1, 40, " MHz");
    adj.drawLabel("t_pulse1 (90°)", params.t_pulse1, 90, " µs");
    adj.drawLabel("t_pulse2 (180°)", params.t_pulse2, 140, " µs");
    adj.drawLabel("T2* (Inhomo)", params.T2_star, 190, " µs");
    adj.drawLabel("T2 (True)", params.T2, 240, " µs");
}

// ============================================
// 物理シミュレーションクラス
// ============================================
class EnsembleSimulator {
    constructor(omega_1, t_pulse1, t_pulse2, T2_star, T2) {
        this.omega_1 = omega_1;
        this.t_pulse1 = t_pulse1;
        this.t_pulse2 = t_pulse2;
        this.T2_star = T2_star;
        this.T2 = T2;
        this.GAMMA = 10;
        this.B_0 = 1;
        this.OMEGA_0 = this.GAMMA * this.B_0; 
        this.B_1 = (PI / 2) / (this.GAMMA * 4.5); 
        this.theta1 = this.GAMMA * this.B_1 * this.t_pulse1;
        this.theta2 = this.GAMMA * this.B_1 * this.t_pulse2;
        this.delta_pulse = this.omega_1 - this.OMEGA_0;
        this.numSpins = 500; 
        this.spins = [];
        let sigma = 1.0 / this.T2_star; 
        for (let i = 0; i < this.numSpins; i++) {
            this.spins.push(randomGaussian(0, sigma));
        }
    }
    calculation_PulseMag() { return 1700 * this.B_1; }
    getSignal(t, isEchoMode, tau) {
        let sumMy = 0;
        for (let i = 0; i < this.numSpins; i++) {
            let dw = this.spins[i]; 
            let total_dw = this.delta_pulse + dw;
            let My = sin(this.theta1);
            let Mx = 0;
            let t_evolve = isEchoMode ? tau : t;
            let phase1 = total_dw * t_evolve;
            let My_1 = My * cos(phase1) + Mx * sin(phase1);
            let Mx_1 = Mx * cos(phase1) - My * sin(phase1);
            if (isEchoMode) {
                // ハードパルス近似：パルス中のオフセットは無視
                let Mz_initial = cos(this.theta1);
                let My_2 = My_1 * cos(this.theta2) - Mz_initial * sin(this.theta2); 
                let Mx_2 = Mx_1; 
                let t_after = t - tau;
                if (t_after > 0) {
                    let phase2 = total_dw * t_after;
                    let My_final = My_2 * cos(phase2) + Mx_2 * sin(phase2);
                    sumMy += My_final;
                }
            } else {
                sumMy += My_1;
            }
        }
        let signal = sumMy / this.numSpins;
        // T2緩和（不可逆）：パルス間隔中のT1回復は無視
        signal *= exp(-t / this.T2);
        return signal * (5000 / PI) * this.B_1;
    }
}

class Canvas {
    constructor() {
        this.canvasWidth = 750; 
        this.canvasHeight = 450;
        this.canvasBackGround = "#000000";
        this.marginLeft = 80;   
        this.marginTop = 10;
        this.marginBottom = 30; 
        this.marginRight = 10;
    }
    create() {
        createCanvas(this.canvasWidth, this.canvasHeight);
        background(this.canvasBackGround);
    }
}

class Adjuster extends Canvas {
    constructor() { super(); }
    
    makeSlider(label, min, max, val, step, yPos) {
        let s = createSlider(min, max, val, step);
        s.position(this.canvasWidth - 220, yPos + 25); 
        s.style('width', '160px'); 
        s.input(updateSimulation); 
        return s;
    }

    drawLabel(label, val, yPos, unit) {
        textAlign(LEFT, TOP);
        textSize(14);
        fill(255);
        text(label, this.canvasWidth - 220, yPos); 
        
        textAlign(RIGHT, TOP);
        fill("#00ffff"); 
        let displayVal = Number.isInteger(val) ? val + ".0" : val.toFixed(1);
        text(displayVal + unit, this.canvasWidth - 60, yPos); 
    }
}

class Graph extends Canvas {
    constructor() {
        super();
        this.axisColor = "#aaaaaa"; 
        this.X_0 = this.marginLeft;
        this.Y_0 = this.canvasHeight / 3 * 2;
        this.GraphWidth = this.canvasWidth - 250; 
    }
    create() {
        stroke(this.axisColor); strokeWeight(1);
        line(this.X_0, this.Y_0, this.GraphWidth, this.Y_0); 
        line(this.X_0, 10, this.X_0, this.canvasHeight - 10); 
        
        fill(200); noStroke(); 
        textAlign(CENTER); textSize(11); 

        for (let i = 0; i < (this.GraphWidth - this.X_0) / 30; i++) {
            let x = this.X_0 + i * 30;
            stroke(100); line(x, this.Y_0 - 5, x, this.Y_0 + 5);
            if (i % 2 == 0) {
                noStroke();
                text(i * 5, x, this.Y_0 + 15);
            }
        }
        textAlign(RIGHT); textSize(12); fill(255);
        text("Time [µs]", this.GraphWidth, this.Y_0 + 35);
        
        textAlign(LEFT); 
        text("Voltage [mV]", 10, 20); 
    }
    
    drawPulse(mag, t, t_start) {
        stroke("#00ff00"); strokeWeight(2); 
        let x = this.X_0 + t_start * (30/5);
        let w = t * (30/5);
        let h = mag * (40/10);
        noFill();
        beginShape();
        vertex(x, this.Y_0);
        vertex(x, this.Y_0 - h);
        vertex(x + w, this.Y_0 - h);
        vertex(x + w, this.Y_0);
        endShape();
    }

    drawFID(sim, t_pulse1, tau) {
        let x_start = this.X_0 + t_pulse1 * (30/5);
        let x_end = this.X_0 + tau * (30/5);
        stroke(255); strokeWeight(1.5);
        let prevX = x_start;
        let prevY = this.Y_0 - sim.getSignal(t_pulse1, false, tau) * 4; 
        for (let x = x_start; x < x_end; x += 2) { 
            let t = (x - this.X_0) / (30/5);
            let y = this.Y_0 - sim.getSignal(t, false, tau) * 4;
            line(prevX, prevY, x, y);
            prevX = x; prevY = y;
        }
    }

    drawSE(sim, t_pulse1, tau) {
        let x_start = this.X_0 + (tau + 2 * t_pulse1) * (30/5);
        let x_end = this.GraphWidth; 
        stroke(255); strokeWeight(1.5);
        let prevX = x_start;
        let t_start_calc = tau + 2 * t_pulse1;
        let prevY = this.Y_0 - sim.getSignal(t_start_calc, true, tau) * 4;
        for (let x = x_start; x < x_end; x += 2) {
            let t = (x - this.X_0) / (30/5);
            let y = this.Y_0 - sim.getSignal(t, true, tau) * 4;
            line(prevX, prevY, x, y);
            prevX = x; prevY = y;
        }
    }
}