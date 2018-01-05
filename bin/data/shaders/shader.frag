#version 150

// デフォルト設定
uniform float u_time;
uniform vec2 u_resolution;
out vec4 outputColor;

const int MAX_MARCHING_STEPS = 255; // ループ回数。この回数分レイを進める
const float MIN_DIST = 0.0; // レイの最短距離 // レイの初期位置
const float MAX_DIST = 100.0; // レイの最大距離
const float EPSILON = 0.0001; // ０に限りなく近い数

const float PI = 3.1415926;
const int oct = 8;
const float per = 0.5;

// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <effect functions>  //


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// function to create random value
// 簡単な方法。raymarchingでもpost-processingでの乱数利用であれば、スクリーンに描画するだけだから、引数はvec2型で、スクリーンの正規化された座標をぶち込めばよい。
float easy_random (vec2 p) {
    return fract(sin(dot(p, vec2(12.9898, 4.1414))) * 43758.5453);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// function of white noise effect
// 単純にpost-processingとしてノイズをスクリーンに走らせる
vec4 white_noise_effect(float rand) {
    return vec4(vec3(rand), 1.0);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// function of night scope effect
vec4 night_scope_effect(vec3 color, vec2 st, vec3 effect_color) {
    
//    vec2 p = st * 2.0 - 1.0;
    
    // colorの平均値
    float dest = (color.r + color.b + color.g) / 3.0;
    
    // 四隅が暗くなるように
    float vignette = sin(u_time) - length(st);
    // stにmodを掛けて複製している場合は、少し変更する必要あり。それでも動作がおかしい。
    // float vignette = 1.0 - length(st*2.0);
    dest *= vignette;
    
    // ノイズをかける
    float noise = easy_random(st+mod(u_time, 10.0));
    dest *= noise * 0.5 + 0.5;
    
    // 走査線を走らせる
    // st.yだから横線が走る。st.y*20で振動数を増やしているから、かける量を小さくすると、線の間隔が大きくなる
    float scanLine = abs(sin(st.y * 50.0 + u_time * 25.0)) * 0.5 + 0.75;
    dest *= scanLine;
    
    // effect_colorはお好みで好きな色を入れる
    return vec4(vec3(dest), 1.0) * vec4(effect_color, 1.0);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// function of value noise effect
float interpolate (float a, float b, float x) {
    float f = (1.0 - cos(x * PI)) * 0.5;
    return a * (1.0 - f) + b * f;
}

float irnd(vec2 p) {
    vec2 i = floor(p);
    vec2 f = fract(p);
    vec4 v = vec4(easy_random(vec2(i.x,     i.y    )),
                  easy_random(vec2(i.x+1.0, i.y    )),
                  easy_random(vec2(i.x,     i.y+1.0)),
                  easy_random(vec2(i.x+1.0, i.y+1.0)));
    return interpolate(interpolate(v.x, v.y, f.x), interpolate(v.z, v.w, f.x), f.y);
}

float noise (vec2 p) {
    float t = 0.0;
    for (int i = 0; i < oct; i++) {
        float freq = pow(2.0, float(i));
        float amp = pow(per, float(oct - i));
        t += irnd(vec2(p.x / freq, p.y / freq)) * amp;
    }
    return t;
}

float snoise(vec2 p, vec2 q, vec2 r) {
    return noise(vec2(p.x,     p.y    )) *     q.x *      q.y +
           noise(vec2(p.x,     p.y+r.y)) *     q.x * (1.0-q.y)+
           noise(vec2(p.x+r.x, p.y    )) *(1.0-q.x)*      q.y +
           noise(vec2(p.x+r.x, p.y+r.y)) *(1.0-q.x)* (1.0-q.y);
}

vec4 value_noise() {
    vec2 p = gl_FragCoord.st + u_time * 100.0;
    float n = snoise(p, gl_FragCoord.st / u_resolution, u_resolution);
    return vec4(vec3(n)*2.50, 1.0);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// 奇数行にだけ引数で指定したcolorのeffectをかける
vec4 odd_Row_Effect(vec3 color) {
    // gl_FragCoordは0.5~resolutionのサイズ（1280とか）を保持する
    // 左下が(0.5, 0.5)だから、-0.5を加えて、0.0にする
    // modは第一要素を第二要素で割った余りだから、2で割った余りが１なら奇数。奇数なら1を返し、それ以外は0.
    // 奇数の時だけeffectを入れる
    bool isOdd = mod(gl_FragCoord.x - 0.5, 2.0) == 1.0;
    if (isOdd) return vec4(color, 1.0);
    else return vec4(1.0);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// pは-1~1の座標値を入れる。
// offsetはずらす値
vec4 createMask(vec2 p, vec2 offset) {
    // aspect比を考慮する場合はコメントアウト
    //    p.y /= u_resolution.x / u_resolution.y;
    
    // vignetteでは1.0からlengthで指定した原点からの距離を引くので、この場合だと、1.0を半径とした円が描かれる。length1.0以上だとマイナスになるから、画面が黒くなる
    // clampでは第一要素を第二、第三引数で指定した値の間に抑え込む。この場合だと0.0~1.0に
    float vignette0 = clamp(1.0 - length(p+offset), 0.0, 1.0); // 右
    float vignette1 = clamp(1.0 - length(p-offset), 0.0, 1.0); // 左
    
    // smoothstep用の値。smoothさせる領域。0.5以下は0に。0.55以上は1.0にする。0.5~0.55の間にかけて滑らかに0から1に移っていく。
    float startD = 0.50;
    float endD = 0.55;
    
    // maskを作成.
    // smoothstepでは第三引数が第一引数から第二引数の間にかけて0.0~1.0に移る。第一引数以下では0、第二引数以上では1.0を返す
    // この場合だとvignetteが0.5~0.55の値になる領域にだけ0.0~1.0に滑らかに移るというeffectをかける
    float value = 0.0;
    value += smoothstep(startD, endD, vignette0);
    value += smoothstep(startD, endD, vignette1);
    
    // clampにすることで二つのマスクが重なっても、明るくなることはない
    return vec4(vec3(clamp(value, 0.0, 1.0)), 1.0);
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //
// 水玉いっぱいeffect
vec4 mizutama (vec2 st) {
    // マイナスで-0.5~0.5に
    vec2 mod_st = mod(st*5.0, 1.0)-0.5;
    
    vec4 color;
    
    // １行おきにエフェクト変えたい場合はここら辺をいじくる
    if (mod(st.x*5.0, 2.0) >= 1.0) {
        float v = clamp(abs(cos(u_time)) - length(mod_st), 0.0, 1.0);
        float vv = smoothstep(0.5, 0.55, v);
        color = vec4(vec3(vv), 1.0);
    } else {
        float v = clamp(abs(sin(u_time)) - length(mod_st), 0.0, 1.0);
        float vv = smoothstep(0.5, 0.55, v);
        color = vec4(vec3(vv), 1.0);
    }

    return color;
}





// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <function of transformation>  //
// thetaにはradians()を通してから代入する

// x軸で回転
mat3 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3 (
                 vec3(1, 0, 0),
                 vec3(0, c, -s),
                 vec3(0, s, c)
                 );
}

// y軸で回転
mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3 (
                 vec3(c, 0, s),
                 vec3(0, 1, 0),
                 vec3(-s, 0, c)
                 );
}

// z軸で回転
mat3 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3 (
                 vec3(c, -s, 0),
                 vec3(s, c, 0),
                 vec3(0, 0, 1)
                 );
}

// xyz同時に回転
// axisはどの軸にどれだけ回転させたいか. もし(1.0, 0.5, 0.0)だとx軸に対して100%y軸に対して50%ということになる
mat3 rotate(float theta, vec3 axis) {
    vec3 a = normalize(axis);
    float c = cos(u_time);
    float s = sin(u_time);
    float r = 1.0 - c;
    return mat3 (
                 a.x * a.x * r + c,
                 a.y * a.x * r + a.z * s,
                 a.z * a.x * r - a.y * s,
                 a.x * a.y * r - a.z * s,
                 a.y * a.y * r + c,
                 a.z * a.y * r + a.x * s,
                 a.x * a.z * r + a.y * s,
                 a.y * a.z * r - a.x * s,
                 a.z * a.z * r + c
                 );
}


// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <function to draw multiple objects at one time> //

// intersetSDFでは、重なった部分のみを描写する
// 前提として、negative値または限りなく０に近い値をdistance functionとして返すと、そこがオブジェクトの表面として画面に描かれる
// 両方がnegative値または限りなく０に近い値を示す場合、それはまさに両方が重なっている場合、のみmax関数でも小さい値を返すことになるので描画される。
float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}

// unionSDFでは字の通り複数オブジェクトを合体させる
// オブジェクトが描かれるのは限りなく値が小さいときというのを思い出すと、なぜmin関数で複数オブジェクトが描かれるのかはわかる
float unionSDF(float distA, float distB) {
    return min(distA, distB);
}

// differenceSDFでは複数オブジェクトの差分を利用して描画する
// max関数だから、両方のオブジェクトが小さい値のみ描画される。
// -distBだから、符合が逆転して、もともとobjectとして描かれた所の値がプラスの大きな値を持って、もともとobjectの外側だった部分の値がマイナスの値を持つことになる。
// その状態で第一引数の値が０に限りなく近い値を返すところでオブジェクトを描くので結果として、第一引数のオブジェクトから第二引数のオブジェクトを引いたようになる
float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

// objectの重なり部分を補間する
// objectの重なり部分はd1,d2とも０付近の値を返す.
// exp(-k * d)の部分ではもしdが0に近いのであれば、1を返すので、1+1で2になる。
// h=2周辺で（詳しく計算していない）いい感じの0付近の数になると予想している。
// -log(h)/kで少しだけマイナス値で引くことでオブジェクトの重なり部分を膨らませている
// あまり深く考えずにobjectの重なり部分だとhが0付近の数になりそれをマイナスで弾くことで重なり部分だけを上手い具合に膨らませることができると思っておけばよい
float smoothMin(float d1, float d2, float k) {
    float h = exp( -k * d1 ) + exp( -k * d2 );
    return -log(h) / k;
}


// ------------------------------------------------------------------------------------------------------------------------------------- //
// sceneSDF2対応のオブジェクト結合関数たち
// .wで距離だけを比べて、値が小さい方のベクトルを返す
vec4 unionSDF2(vec4 d1, vec4 d2) {
    return d1.w < d2.w ? d1 : d2;
}

// .wで距離だけを比べて、値が大きい方のベクトルを返す
vec4 intersectSDF2(vec4 d1, vec4 d2) {
    return d1.w > d2.w ? d1 : d2;
}

// .wで距離だけを比べて、値が大きい方のベクトルを返す
vec4 differenceSDF2(vec4 d1, vec4 d2) {
    return d1.w > -d2.w ? d1 : vec4(d2.rgb, -d2.w);
}

// 距離の部分だけsmoothMinを適用。何か他にもやり方あるような気もする
vec4 unionWithSmoothSDF(vec4 d1, vec4 d2, float k) {
    float smooth_d = smoothMin(d1.w, d2.w, k);
    vec4 s = d1.w < d2.w ? d1 : d2;
    s.w = smooth_d;
    return s;
}


// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <object function> //


// simple sphere
// この関数が0.0を返せばsphereの表面。+ならoutside,-ならinside。
// 1(sphereSize) = x^2 + y^2 + z^2 を x^2 + y^2 + z^2がlength()でそこから1を右辺に移動させ、=0で方程式が成り立つ場合にsphereを描く
// sphereの中心を原点として、そこから半径1(sphereSize)をとる球体を置いて、各レイの先端位置との距離を測る
float sphereSDF(vec3 samplePoint, float sphereSize) {
    return length(samplePoint) - sphereSize;
}


// distance function of cube
float cubeSDF(vec3 p, vec3 cubeSize) {
    // 各点の絶対値からcubeの大きさをひく（cubeの大きさは実際にはcubeSizeの二倍になる）
    // insideは0またはマイナス値なる。outsideは0またはプラス値になる
    vec3 d = abs(p) - cubeSize;
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));
    return insideDistance + outsideDistance;
}

// round cube
// 基本上のcubeと同じで、角を丸めるために変数roundを用意して、最終結果でround分を引く
float cubeSDF_round(vec3 p, vec3 cubeSize) {
    // round : どれだけ丸みを帯びさせるか // この値だけ外側に膨らませる(その分オブジェクトが大きくなる)
    float round = 0.1;
    vec3 d = abs(p) - cubeSize;
    float insideDistance = min(max(d.x, max(d.y, d.z)), 0.0);
    float outsideDistance = length(max(d, 0.0));
    return (insideDistance + outsideDistance) - round;
}

// distance function of cylinder
// h : height
// r : radius
float cylinderSDF( vec3 p, float h, float r) {
    // 単純にxy座標において、半径rの円内にあるかどうか。半径ないならマイナス値で半径外ならプラス値になる。円上なら0になる
    float inOutRadius = length(p.xy) - r;
    // 単純にz軸の絶対値をとることで、原点からの純粋な距離を測りそこから円柱のHeightの半分だけ引く。(半分なのは原点を挟んでいて、z軸の絶対値から引いているから)
    float inOutHeight = abs(p.z) - h/2.0;
    // もし円柱の内側なら両方がマイナス値なので、そのマイナス値が返される。他は０になる
    float insideDistance = min(max(inOutRadius, inOutHeight), 0.0);
    // もし円柱の内側なら両方がマイナス値なので、max関数で０になる。もしoutsideなら０以上のプラス値が返させる。
    float outsideDistance = length(max(vec2(inOutRadius, inOutHeight), 0.0));
    // 0もしくは,限りなく０に近い値を返す点において、オブジェクトを描く
    return insideDistance + outsideDistance;
}

// distance function of eacy clinder
// 両端が切れないcylinder
// 操作性低い
float easyCylinderSDF (vec3 p) {
    // 0にセットすると原点からの距離を測れる
    vec2 c = vec2(0.0, 0.0);
    // radius
    float radius = 0.5;
    // ここでpのxyzの中の２つだけしか使わないことで両端が切れないcylinderを描ける
    // 単純に原点からの距離を測ってradiusで引く。その値が限りなく0に近い値であれば、そこが表面
    return length(p.yz - c.xy) - radius;
}

// distance function of the plane
// 内積を使う。
// 例えば画面下に床を描きたいのであれば、nはvec4(0.0, 1.0, 0.0, 1.0)にする
// 内積でx, zは無視されy軸に対して、レイのy座標がマイナスであれば、内積結果はマイナスになる.
// それでn.wを足し合わせた結果が0付近になればその集合が描画されるところである。
float planeSDF(vec3 p, vec4 n) {
    // n must be normalized
    vec3 nn = normalize(n.xyz);
    return dot(p, nn) + n.w;
}


// distance function of the torus
// t.xはtorusの中心からどれだけ距離をおいてパイプを作るか
// t.yはパイプの太さ。この値が大きいとパイプは太くなる。
// p.xyならz軸方向を向く。p.yzならx軸方向を向く。p.xzならy軸方向を向く
// rの第一要素では単純に、この場合だとxy平面で円を描くのと同じで円の円周で０を返すようにしている。
// rで円周上であるならrの第一要素は0だから、最後のreturnではp.zの距離がt.y(パイプの太さ)と同じところで０を返すようになる
// 変更：引数 vec3 p -> vec2 base, float up に変更。こうすることで、向きの異なるtorusを作れる。 baseはどの平面をベースとするか。upはどちらを上としてtorusを描くか。やっていることは前と同じで単にバラしただけ。
float torusSDF(vec2 base, float up, vec2 t ) {
    vec2 r = vec2( length( base ) - t.x, up );
    return length( r ) - t.y;
}

// distance function of the cone
// 少しばかり難しい。
// returnのところは円柱の描画と全く同じ
// d1はconeの高さを測っている。
// d2はmax関数なので、q.yがマイナスつまりy<0で、dotの計算が大きくなれば良い。
// dotが０に近い値を返すのは、二つのベクトルが90度の角度で離れている時。
// だからc.xを大きくすれば円錐の円は大きくなり、c.yを大きくすれば円錐は細長くなる。
float coneSDF( vec3 p, vec3 c ) {
    vec2 q = vec2( length(p.xz), p.y );
    float d1 = -q.y - c.z;
    float d2 = max( dot(q,c.xy), q.y );
    float insideDistance = min( max(d1,d2), 0.0 );
    float outsideDistance = length( max(vec2(d1,d2),0.0) );
    return insideDistance + outsideDistance;
}

// distance function of the circle
// 単に円柱を押しつぶしただけだから、円柱と全く同じアルゴリズム。
float circleSDF( vec3 p, float d ) {
    float d1 = length(p.xy) - d;
    float d2 = abs(p.z) - .001;
    float inside = min(max(d1, d2), 0.0);
    float outside = length(max(vec2(d1, d2), 0.0));
    return inside + outside;
}

// distance function of the rectangle
// 単に四角形を描く。これもcubeを潰しただけだから、アルゴリズムはcubeと同じ。
float rectangleSDF( vec3 p, vec2 rectSize) {
    vec2 d1 = abs(p.xy) - rectSize;
    float d2 = abs(p.z) - 0.001;
    float inside = min(max(d2, max(d1.x, d1.y)), 0.0);
    float outside = length(max(vec3(d1, d2), 0.0));
    return inside + outside;
}


// この関数をdistance_functionのハブとしておくことで、複数オブジェクトを描いたり、重ねて描いたりすることを容易にする。
float sceneSDF(vec3 samplePoint) {
    float cylinderRadius = 0.3;
    float cylinder1 = cylinderSDF( samplePoint, .50, cylinderRadius );
    float cylinder2 = cylinderSDF( rotateX(radians(90.0)) * samplePoint, 2.0, cylinderRadius );
    float cylinder = smoothMin (cylinder1, cylinder2, 7.0);
    
    float floor = planeSDF(samplePoint, vec4(0.0, 0.0, 1.0, .10));
    return unionSDF(cylinder, floor);
}

// 各オブジェクトが固有のカラーを持てるように変更
// returnをvec4にしている。んで各オブジェクトにカラー.rgbでセット
vec4 sceneSDF2(vec3 samplePoint) {
    
    vec3 pos_sp = samplePoint;
    float sp = sphereSDF(pos_sp, 1.0);
    vec4 sphere = vec4(vec3(0.0), sp);
    
    vec3 pos_c = samplePoint;
    float c = cylinderSDF(pos_c, 2.5, 0.5);
    vec4 cylinder = vec4(vec3(1.0, 0., 0.), c);
    
    
    
    return differenceSDF2(sphere, cylinder);
}



// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <basic function> //

// ここでレイを作る。
// xy方向は単に -u_resolution/2 ~ u_resolution/2 に座標変換したスクリーンのxyを入れる。
// z方向に関してはfieldOfViewの角度に応じて、z軸方向へのレイの大きさを変える
vec3 rayDirection(float fieldOfView) {
    vec2 xy = gl_FragCoord.xy - u_resolution / 2.0;
    float z = u_resolution.y / tan(radians(fieldOfView)/2.0);
    return normalize(vec3(xy, -z));
}

vec3 rayDirection2(vec2 xy, float fieldOfView) {
    float fov = fieldOfView * (PI / 180.0);
    return normalize(vec3(sin(fov)*xy.x, sin(fov)*xy.y, -cos(fov)));
}

// ここでレイを実際に飛ばしてオブジェクトまでの距離を測る。
// eye : eye point。レイの初期位置にeyeと仮定しておく。(人間にとってその方がわかりやすいだけ)
// marchingDirection : レイの方向。ここは単にレイを入れるだけ。
// start : 初期状態でレイがカメラからどれだけ離れているか
// end : カメラからのレイの最大距離。これ以上のレイの値を追うことはしない。
float shortestDistanceToSurface(vec3 eye, vec3 marchingDirection, float start, float end) {
    float depth = start;
    for ( int i = 0; i < MAX_MARCHING_STEPS; i++ ) {
        float rayPos = sceneSDF2( eye + depth * marchingDirection ).w;
        // EPSILONより小さい、つまりオブジェクトの表面だとわかり次第をreturnでdepthを返してループを抜ける
        // rayPosとdepthは別物。depthは次第に大きくなるけど、rayPosはオブジェクトとの交点に近くなるにつれて小さくなる（だって、オブジェクトとの距離を示す変数だから）
        if ( rayPos < EPSILON ) {
            return depth;
        }
        // レイを進める
        depth += rayPos;
        // depthが最大距離を超えるとendを返してループを抜ける
        if ( depth >= end ) {
            return end;
        }
    }
    // returnされなかったものに対してもendを返す
    return end;
}


// floatではなくベクトルで返す
// ちょっぴり面倒なのは、endもベクトルで返さないといけないこと
vec2 shortestDistanceToSurface2(vec3 eye, vec3 marchingDirection, float start, float end) {
    // depthをベクトルとして扱う。
    // 第一要素にはレイの深さ。つまりこの値はどんどん大きくなる
    // 第二要素はレイの先端位置（オブジェクトとの距離）。objectと交差しないのであれば値は開いていくが、もし交差するのであればどんどんobjectに近づくわけだから、値は小さくなる
    vec2 depth;
    vec2 max = vec2(end);
    depth.x = start;
    for ( int i = 0; i < MAX_MARCHING_STEPS; i++ ) {
        depth.y = sceneSDF2( eye + depth.x * marchingDirection ).w;
        
        if ( depth.y < EPSILON ) {
            // ベクトルdepthにまとめることで、レイの深さとobjectとの距離の両方を返せる。
            // このdepthはEPSILON以下であることを条件として返すわけだから、非常に小さい(極めて０に近い)値を返す
            return depth;
        }
        // レイを進める
        depth.x += depth.y;
        // depthが最大距離を超えるとendを返してループを抜ける
        if ( depth.x >= max.x ) {
            return max;
        }
    }
    // returnされなかったものに対してもendを返す
    return max;
}

// SDFの勾配を求めて、各ポイントにおける法線を算出。
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
                          sceneSDF2(vec3(p.x + EPSILON, p.y, p.z)).w - sceneSDF2(vec3(p.x - EPSILON, p.y, p.z)).w,
                          sceneSDF2(vec3(p.x, p.y + EPSILON, p.z)).w - sceneSDF2(vec3(p.x, p.y - EPSILON, p.z)).w,
                          sceneSDF2(vec3(p.x, p.y, p.z + EPSILON)).w - sceneSDF2(vec3(p.x, p.y, p.z - EPSILON)).w
    ));
}

// the function for generating shadow
//
float genShadow(vec3 refpos, vec3 raydir) {
    float h = 0.0;
    float c = 0.01;
    float r = 1.0;
    float shadowCoef = 0.5;
    for (float t = 0.0; t < 50.0; t++) {
        h = sceneSDF2( refpos + raydir * c ).w;
        if ( h < 0.001 ) {
            return shadowCoef;
        }
        r = min(r, h * 16.0 / c);
        c += h;
    }
    return 1.0 - shadowCoef + r * shadowCoef;
}

// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <lighting function> //

// k_d : diffuse color
// k_s : specular color
// alpha : shininess coefficient
// p : position of point begin lit
// eye : position of the camera
// lightPos : the position of the light
// lightIntensity : color/light intensity of the lihgt
vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye, vec3 lightPos, vec3 lightIntensity) {
    vec3 N = estimateNormal(p); // N : Normal
    vec3 L = normalize(lightPos-p); // L : pから光源方向へのベクトル
    vec3 R = normalize(reflect(-L, N)); // R : 反射ベクトル（光源から点pに向かって放たられる光に対する反射）
    vec3 V = normalize(eye-p); // pから目線（カメラ位置）方向へのベクトル
    
    float dotLN = dot(L, N); // ベクトルLとNの内積を計算
    float dotRV = dot(R, V); // ベクトルRとVの内積を計算
    
    if ( dotLN < 0.0 ) {
        // もし内積が０以下、つまり二つのベクトルが９０以上開いていたらライトを消す（0を返す）
        return vec3 (0.0, 0.0, 0.0);
    }
    if ( dotRV < 0.0 ) {
        // pから目線方向へのベクトルと反射ベクトルの角度が９０以上開いていたらdiffuseのみを適用する
        return lightIntensity * (k_d*dotLN);
    }
    return lightIntensity * (k_d*dotLN+k_s*pow(dotRV, alpha));
}

// vec3で返される値は反射後のRGBカラーの値。
// k_a : ambient color
// k_d : diffuse color
// k_s : specular color
// alpha : shininess coefficient. この定数が大きいほど、鏡面ハイライトが小さく強くなる
// p : position of point beging lit
// eye : position of camera
vec3 phongillumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye, vec2 dist ) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    vec3 light1Pos = vec3(4.0*sin(u_time),
                          2.0,
                          4.0*cos(u_time));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    color += phongContribForLight(k_d, k_s, alpha, p, eye, light1Pos, light1Intensity);
    
//    vec3 light2Pos = vec3(2.0 * sin(0.37 * u_time),
//                          2.0 * cos(0.37 * u_time),
//                          2.0);
//    vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
//    color += phongContribForLight(k_d, k_s, alpha, p, eye, light2Pos, light2Intensity);
    
    
    float shadow = 1.0;
    //-------------------------------generating tile-------------------------------------//
    //------------------------------generating shadow------------------------------------//
    
//    if ( abs(dist.y) < EPSILON ) {
//        vec3 normal = estimateNormal(p);
//        float diff = clamp(dot(light1Pos, normal), 0.5, 1.0);
//
//        shadow = genShadow(p + normal*0.001, light1Pos);
//
//        float u = 1.0 - floor(mod(p.x, 2.0));
//        float v = 1.0 - floor(mod(p.z, 2.0));
//        if (( u == 1.0 && v < 1.0 ) || ( u < 1.0 && v == 1.0 )) {
//            diff *= .5;
//        }
//        color = vec3(1.0)*diff;
//    }
    
    //----------------------------------------------------------------------------------//
    
    return color * max(shadow, 0.5);
}




// viewMatrixを作る。カメラ中心の座標にする。
mat4 viewMatrix(vec3 eye, vec3 center, vec3 up) {
    vec3 f = normalize( center - eye );
    vec3 s = normalize( cross(f, up) );
    vec3 u = cross(s, f);
    return mat4 (
        vec4(s, 0.0),
        vec4(u, 0.0),
        vec4(-f, 0.0),
        vec4(0.0, 0.0, 0.0, 1.0)
    );
}






// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// ------------------------------------------------------------------------------------------------------------------------------------- //
// <main> //

void main () {
    // スクリーン座標を上下左右を-1~1にする (左下が-1, -1で右上が1,1)
    vec2 st = (gl_FragCoord.xy * 2.0 - u_resolution.xy) / min(u_resolution.x, u_resolution.y);
    
    // modで複製するなら
//    st = mod(st*2.0, 1.0)-0.5;
    
    // make uv coordinate
    // -1~1のuv座標を作る
//    vec2 uv = (gl_FragCoord.xy * 2.0 - u_resolution.xy) / min(u_resolution.x, u_resolution.y);
//    if (uv.y != 0.0) {
//        uv.x /= abs(uv.y*2.0);
//    }
//    uv.y = abs(uv.y);
////    uv.x += u_time;
//    uv = mod(uv*5.0, 1.0)-0.5;
    
    
    // fieldOfViewの角度を渡してレイを作成
    vec3 viewDir = rayDirection(45.0);
//    vec3 viewDir = rayDirection2(st, 45.0);
    
    // カメラの位置を決める
//    vec3 eye = vec3(8.0, sin(u_time*0.2)*5.0, 7.0);
//    vec3 eye = vec3(8.0, 4.0, 5.0);
    vec3 eye = vec3(cos(u_time)*10.0, 4.0, sin(u_time)*10.0);
//    vec3 eye = vec3(0.0, 0.0, 10.0);
    
    // viewMatrixを作る。ここでカメラ中心の座標にする(openGLの行列チュートリアル見ると分かりやすい)
    mat4 viewToWorld = viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
    // ここで上でつくったカメラ中心の座標にたいして、各ピクセルに放たれるレイのベクトルを掛け合わせる。
    // 要するにここら辺ではカメラに対応するように座標変換している
    vec3 worldDir = (viewToWorld * vec4(viewDir, 0.0)).xyz;
    
    vec2 dist = shortestDistanceToSurface2( eye, worldDir, MIN_DIST, MAX_DIST );
    
    
    // shortestDistanceToSurface関数でレイを進めていることと同じ。
    // 違いはvec3でそのオブジェクトの表面のベクトル情報を取っていること。
    // surfPosは単にオブジェクトとの交差点の座標を持つベクトルである。
    vec3 surfPos = eye + dist.x * worldDir;
    
    vec3 K_a = sceneSDF2(surfPos).rgb;
    vec3 K_d = sceneSDF2(surfPos).rgb;
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;
    
    vec3 color = phongillumination( K_a, K_d, K_s, shininess, surfPos, eye, dist );
    
    
    // distではoutsideならMAX_DISTが入る。(shortestDistanceToSurfaceでMAX_DISTが返されている)
    // だからMAX_DISTから小さな値を引いた数より大きい場合は、ピクセルを黒で塗りつぶす。
    // dist >= MAX_DISTでも同じ効果が得られる。
    if( dist.x > MAX_DIST - EPSILON ) {
        outputColor = vec4(1.0);
        
        // post-effectをかけたい時、スクリーン全体にpost-effectをかけたい。
        // そうするために、レイが衝突していない領域でもpost-effectをかける必要あり。
        // 衝突していない領域をvec4(0.0)とするなら、addする。
        // 衝突していない領域をvec4(1.0)とするなら、multiplyする。
//        outputColor *= night_scope_effect(color, st, vec3(0.9, 0.9, 0.7));
//        outputColor *= white_noise_effect(easy_random(st+u_time));
//        outputColor *= odd_Row_Effect(vec3(1., 0., 0.));
        return;
    }
    
    outputColor = vec4(color, 1.0);
    
    // from here you can try some post-effects
//    outputColor *= night_scope_effect(color, st, vec3(0.9, 0.9, 0.7));
//    outputColor *= odd_Row_Effect(vec3(1., 0., 0.));
//    outputColor *= mizutama(st);
    
}
