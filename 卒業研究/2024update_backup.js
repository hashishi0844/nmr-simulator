/** JavaScript で使われる関数などの説明 */
// PI...円周率
// sqrt()...平方根
// pow(底, 指数)...累乗
// randomGaussian(mean, sd)...正規分布に従う乱数を生成

/** グローバル変数（シミュレーション用） */
// 初期値だけここで設定し、あとはスライダーが管理します
let params = {
    omega_1: 10.0,
    t_pulse1: 4.5,
    t_pulse2: 9.0,
    T2_star: 5.0,
    T2: 100.0
};
const tau = 30;

// スライダーの入れ物
let sliders = {};

// シミュレーターのインスタンス
let mySim;

/** ページ読み込み時の処理（最初に1回だけ実行） */
function setup() {
    new Canvas().create();
    
    // スライダー等のUIを作成（1回だけ）
    createInterface();

    // 最初の描画を実行
    updateSimulation();
}

/** 画面の更新処理（スライダーを動かすたびに呼ばれる） */
function updateSimulation() {
    // 1. 背景をリセット（真っ黒に塗りつぶす）
    background("#000000");

    // 2. スライダーから現在の値を読み取る
    params.omega_1  = sliders.omega_1.value();
    params.t_pulse1 = sliders.t_pulse1.value();
    params.t_pulse2 = sliders.t_pulse2.value();
    params.T2_star  = sliders.T2_star.value();
    params.T2       = sliders.T2.value();

    // 3. 値の数値を画面に表示する（ラベル描画）
    drawLabels();

    // 4. シミュレーターを計算（新しい値でアンサンブル生成）
    mySim = new EnsembleSimulator(
        params.omega_1, 
        params.t_pulse1, 
        params.t_pulse2, 
        params.T2_star, 
        params.T2
    );
    
    // 5. グラフとパルスの描画
    let graph = new Graph();
    graph.create(); // 軸を描く

    const pulseMag = mySim.calculation_PulseMag();
    graph.drawPulse(pulseMag, params.t_pulse1, 0);
    graph.drawPulse(pulseMag, params.t_pulse2, tau);
    
    graph.drawFID(mySim, params.t_pulse1, tau);
    graph.drawSE(mySim, params.t_pulse1, tau);
}


/** 自動的に60回/sで繰り返し処理される */
function draw() {  
    noLoop(); // ループしない（スライダー操作時のみ更新）
}


// ============================================
// UI作成・管理クラス
// ============================================
function createInterface() {
    let adj = new Adjuster();
    
    // createSlider(最小, 最大, 初期値, ステップ)
    sliders.omega_1  = adj.makeSlider("omega_1", 7.5, 12.5, params.omega_1, 0.1, 40);
    sliders.t_pulse1 = adj.makeSlider("t_pulse1", 0.0, 10.0, params.t_pulse1, 0.1, 90);
    sliders.t_pulse2 = adj.makeSlider("t_pulse2", 0.0, 20.0, params.t_pulse2, 0.1, 140);
    sliders.T2_star  = adj.makeSlider("T2* (Inhomo)", 1.0, 20.0, params.T2_star, 1.0, 190);
    sliders.T2       = adj.makeSlider("T2 (True)", 10.0, 500.0, params.T2, 10.0, 240);
}

function drawLabels() {
    // スライダーの横や上に数値を表示する
    let adj = new Adjuster();
    adj.drawLabel("omega_1", params.omega_1, 40);
    adj.drawLabel("t_pulse1 (90°)", params.t_pulse1, 90);
    adj.drawLabel("t_pulse2 (180°)", params.t_pulse2, 140);
    adj.drawLabel("T2* (Inhomo)", params.T2_star, 190);
    adj.drawLabel("T2 (True)", params.T2, 240);
}


// ============================================
// 物理シミュレーションクラス (Ensemble Model)
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
            let inhomogeneity = randomGaussian(0, sigma);
            this.spins.push(inhomogeneity);
        }
    }

    calculation_PulseMag() {
        return 1700 * this.B_1;
    }

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
        signal *= exp(-t / this.T2);
        
        return signal * (5000 / PI) * this.B_1;
    }
}


class Canvas {
    constructor() {
        this.canvasWidth = 640;
        this.canvasHeight = 450;
        this.canvasBackGround = "#000000";
        this.marginLeft = 60;
        this.marginTop = 10;
        this.marginBottom = 10;
        this.marginRight = 10;
    }
    create() {
        createCanvas(this.canvasWidth, this.canvasHeight);
        background(this.canvasBackGround);
    }
}

class Adjuster extends Canvas {
    constructor() { super(); }
    
    // スライダーを作る関数
    makeSlider(label, min, max, val, step, yPos) {
        let s = createSlider(min, max, val, step);
        // 位置調整: 画面の右側 (幅640 - 180 = 460くらい)
        s.position(this.canvasWidth - 180, yPos + 25); 
        s.style('width', '140px'); // スライダーの幅
        s.input(updateSimulation); // 動かしたら画面更新
        return s;
    }

    // 数値とラベルを描画する関数
    drawLabel(label, val, yPos) {
        textAlign(LEFT, TOP);
        textSize(14);
        fill(255);
        text(label, this.canvasWidth - 180, yPos); // ラベル
        
        textAlign(RIGHT, TOP);
        fill("#00ffff"); // 数値だけ色を変えて目立たせる
        // 小数点以下の桁数を調整
        let displayVal = val;
        if (Number.isInteger(val)) {
            displayVal = val + ".0"; 
        } else {
            displayVal = val.toFixed(1);
        }
        text(displayVal, this.canvasWidth - 40, yPos); // 数値
    }
}

class Graph extends Canvas {
    constructor() {
        super();
        this.axisColor = "#999999";
        this.X_0 = this.marginLeft;
        this.Y_0 = this.canvasHeight / 3 * 2;
    }
    create() {
        stroke(this.axisColor); strokeWeight(1);
        line(this.X_0, this.Y_0, this.canvasWidth - 10, this.Y_0); 
        line(this.X_0, 10, this.X_0, this.canvasHeight - 10); 
        
        fill(150); noStroke(); textAlign(CENTER); textSize(10);
        for (let i = 0; i < this.canvasWidth / 30; i++) {
            let x = this.X_0 + i * 30;
            if (x > this.canvasWidth) break;
            stroke(100); line(x, this.Y_0 - 5, x, this.Y_0 + 5);
            if (i % 2 == 0) {
                noStroke();
                text(i * 5, x, this.Y_0 + 20);
            }
        }
        textAlign(LEFT);
        text("Time [µs]", this.canvasWidth - 60, this.Y_0 + 35);
        
        textAlign(LEFT);
        text("Voltage [mV]", 10, 20);
    }
    
    drawPulse(mag, t, t_start) {
        stroke("#00ff00"); strokeWeight(2); // パルスは緑色で見やすく
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
        let x_end = this.canvasWidth - 10;
        
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