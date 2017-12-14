#version 150

// デフォルト設定
uniform float u_time;
uniform vec2 u_resolution;
out vec4 outputColor;

const int MAX_MARCHING_STEPS = 255; // ループ回数。この回数分レイを進める
const float MIN_DIST = 0.0; // レイの最短距離
const float MAX_DIST = 100.0; // レイの最大距離
const float EPSILON = 0.0001; // ０に限りなく近い数


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
    return d1.w > -d2.w ? d1 : d2;
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
// object function //


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
    return dot(p, n.xyz) + n.w;
}


// distance function of the torus
// t.xはtorusの中心からどれだけ距離をおいてパイプを作るか
// t.yはパイプの太さ。この値が大きいとパイプは太くなる。
// p.xyならz軸方向を向く。p.yzならx軸方向を向く。p.xzならy軸方向を向く
// rの第一要素では単純に、この場合だとxy平面で円を描くのと同じで円の円周で０を返すようにしている。
// rで円周上であるならrの第一要素は0だから、最後のreturnではp.zの距離がt.y(パイプの太さ)と同じところで０を返すようになる
float torusSDF(vec3 p, vec2 t) {
    vec2 r = vec2(length(p.xy)-t.x, p.z);
    return length(r) - t.y;
}




// この関数をdistance_functionのハブとしておくことで、複数オブジェクトを描いたり、重ねて描いたりすることを容易にする。
float sceneSDF(vec3 samplePoint) {
    float cylinderRadius = 0.3;
    float cylinder1 = cylinderSDF( samplePoint, 2.0, cylinderRadius );
    float cylinder2 = cylinderSDF( rotateX(radians(90.0)) * samplePoint, 2.0, cylinderRadius );
    float cylinder = smoothMin (cylinder1, cylinder2, 7.0);
    
    float floor = planeSDF(samplePoint, vec4(0.0, 0.0, 1.0, .10));
    
    return unionSDF(cylinder, floor);
//    return unionSDF(sphere, box);
}

// 各オブジェクトが固有のカラーを持てるように変更
// returnをvec4にしている。んで各オブジェクトにカラー.rgbでセット
vec4 sceneSDF2(vec3 samplePoint) {
    float torus = torusSDF(samplePoint, vec2(0.75, 0.25));
    vec4 t = vec4(vec3(1.0, 1.0, 1.0), torus);
    return t;
}







// ここでレイを作る。
// xy方向は単に -u_resolution/2 ~ u_resolution/2 に座標変換したスクリーンのxyを入れる。
// z方向に関してはfieldOfViewの角度に応じて、z軸方向へのレイの大きさを変える
vec3 rayDirection(float fieldOfView) {
    vec2 xy = gl_FragCoord.xy - u_resolution / 2.0;
    float z = u_resolution.y / tan(radians(fieldOfView)/2.0);
    return normalize(vec3(xy, -z));
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



// SDFの勾配を求めて、各ポイントにおける法線を算出。
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
                          sceneSDF2(vec3(p.x + EPSILON, p.y, p.z)).w - sceneSDF2(vec3(p.x - EPSILON, p.y, p.z)).w,
                          sceneSDF2(vec3(p.x, p.y + EPSILON, p.z)).w - sceneSDF2(vec3(p.x, p.y - EPSILON, p.z)).w,
                          sceneSDF2(vec3(p.x, p.y, p.z + EPSILON)).w - sceneSDF2(vec3(p.x, p.y, p.z - EPSILON)).w
    ));
}

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
vec3 phongillumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 eye ) {
    const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
    vec3 color = ambientLight * k_a;
    vec3 light1Pos = vec3(4.0*sin(u_time),
                          2.0,
                          4.0*cos(u_time));
    vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
    color += phongContribForLight(k_d, k_s, alpha, p, eye, light1Pos, light1Intensity);
    
    vec3 light2Pos = vec3(2.0 * sin(0.37 * u_time),
                          2.0 * cos(0.37 * u_time),
                          2.0);
    vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
    color += phongContribForLight(k_d, k_s, alpha, p, eye, light2Pos, light2Intensity);
    return color;
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


void main () {
    // スクリーン座標を上下左右を-1~1にする (左下が-1, -1で右上が1,1)
    vec2 st = (gl_FragCoord.xy * 2.0 - u_resolution.xy) / min(u_resolution.x, u_resolution.y);
    
    // fieldOfViewの角度を渡してレイを作成
    vec3 viewDir = rayDirection(45.0);
    // カメラの位置を決める
//    vec3 eye = vec3(8.0, sin(u_time*0.2)*5.0, 7.0);
//    vec3 eye = vec3(8.0, 5.0, 7.0);
    vec3 eye = vec3(0.0, 0.0, 15.0);
    
    // viewMatrixを作る。ここでカメラ中心の座標にする(openGLの行列チュートリアル見ると分かりやすい)
    mat4 viewToWorld = viewMatrix(eye, vec3(0.0, 0.0, 0.0), vec3(0.0, 1.0, 0.0));
    // ここで上でつくったカメラ中心の座標にたいして、各ピクセルに放たれるレイのベクトルを掛け合わせる。
    // 要するにここら辺ではカメラに対応するように座標変換している
    vec3 worldDir = (viewToWorld * vec4(viewDir, 0.0)).xyz;
    
    float dist = shortestDistanceToSurface( eye, worldDir, MIN_DIST, MAX_DIST );
    
    // distではoutsideならMAX_DISTが入る。(shortestDistanceToSurfaceでMAX_DISTが返されている)
    // だからMAX_DISTから小さな値を引いた数より大きい場合は、ピクセルを黒で塗りつぶす。
    // dist >= MAX_DISTでも同じ効果が得られる。
    if( dist > MAX_DIST - EPSILON ) {
        outputColor = vec4(0.0, 0.0, 0.0, 0.0);
        return;
    }
    
    // shortestDistanceToSurface関数でレイを進めていることと同じ。
    // 違いはvec3でそのオブジェクトの表面のベクトル情報を取っていること。
    vec3 surfPos = eye + dist * worldDir;
    
    vec3 K_a = sceneSDF2(surfPos).rgb;
    vec3 K_d = sceneSDF2(surfPos).rgb;
    vec3 K_s = vec3(1.0, 1.0, 1.0);
    float shininess = 10.0;
    
    vec3 color = phongillumination( K_a, K_d, K_s, shininess, surfPos, eye );
    
    outputColor = vec4(color, 1.0);
}
