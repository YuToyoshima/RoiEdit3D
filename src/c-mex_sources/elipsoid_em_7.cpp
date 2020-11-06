/* 更新履歴
elipsoid_em:
ver1: 
    elipsoid_mex_21.cppより派生。
    EM-like optimizationの実装。
    第二引数の扱いを[a,b,c,t,u,v,xc,yc,zc,intensity]'から
    [S_11,S_12,S_13,S_22,S_23,S_33,xc,yc,zc,pi]' へ変更した。
ver2: 
    m1step, m2stepで、各voxelごとに、各楕円体の計算範囲か判定していたが、
    元画像の格納形式をstraight forwardにすることで判定を簡略化した。
    Hstepについては判定が不要になり高速化できた。
    またSTL::vectorを用いて可変長のテーブルを作成し、計算が必要な楕円体の番号を
    各voxelごとに保持するようにして、判定にかかる時間を除いた。
ver3:
    m1step, m2stepで、z方向に飛び飛びのとき、間を飛ばして必要なz位置でだけ計算
    するよう、計算が必要な楕円体の番号を各zごとに保持するテーブルを追加した。
ver4:
    openmp対応
ver5:
    sigmaの対角成分が負値になる場合、sigmaの更新を破棄
ver6:
    sigmaの行列式が負値になる場合、sigmaの更新を破棄
ver7:
    閾値以下のボクセルについて、輝度値を0とみなして残差を計算するように変更。
    



elipsoid_cuda:
ver1: 初期型, 100elipsoid*100parallel*(x=256,y=512,z=153)で50sec強
ver2: reduction kernelをまとめた
ver3: 積極的な枝刈りを導入
ver4: 枝狩り時のinstructionを改良 32sec
ver5: pinned memory 導入
ver6: スレッドブロック数の過多(過剰分割)によるパフォーマンスへの影響を考慮し、
      複数のz座標について同一thread内で扱う方針に変更。31sec
ver7: 最初から必要なpixel数分だけスレッドを起動させることにする。
	  -> ベース座標をdev_param_cuda_intとしてkernelに与えようとしたところ、
	  なぜかdev_param_cuda (float)が先頭から上書きされる事象が発生したため、
	  dev_param_cuda_intを使わずに済む方法へ変更した。
	  -> Z下限から順にZ増分方向へ計算したところ、Z下限でrealmin('single')を
	  下回る値が発生するなどしてNaNが生じたため中断。
ver8: Z方向について、中央付近をbaseとして、Z上下方向へ計算していくことにする。
	base点のintensityを記憶していると桁あふれやゼロへの丸めが発生するので、
	base点のdistanceを代わりに記憶しておくことにする。
	あわせて、scaleをlog(scale)としておいて、distへ取り込んでdist_scaleとすれば
	各z点における計算量は増えない。
	intensity = scale * exp(-dist)
			  = exp(log(scale)-dist)
			  = exp(dist_scale)
ver9: 楕円体の輝度が1を超えないもの(dist_scale<0)についてはexpf以下を省く。
	cudaMemcpyを1回で終わらせるためim_orig (pinned memory)を画像8個分確保する。
	myreduce関数を1回で終わらせるため改変。 
    100elipsoid*100parallel*(x=256,y=512,z=153)で2.2sec程度になった。
    myreduce関数の起動数を最適化したら1.8sec程度になった。
ver10: openmp対応 
    -> 速度はあまり変わらず。
ver11: C++側のベンチマークの結果、cudaMemCpy(im_orig)、elipsoid, myreduce2に
    それぞれ 110msec, 4msec, 30msecかかること、im_origのサイズに応じて
    cudaMemcpyの所要時間がスケールすることがわかった。
    -> im_origの転送時間を削るため、im_orig1つ分を初期化時にdev側へ転送しておき、
    計算時にkernel関数を使って各im_diffへ展開すれば高速化できるはず。
    cudaMemcpyをmemcpyD2Dにしたところ110msec->8msecへ短縮できた。
    100elipsoid*100parallel*(x=256,y=512,z=153)で500msec程度になった。
    -> memcpyD2Dとmyreduce2は内容が似ていて、cuda profilerではそれぞれ
    8msecと6msec程度なのに、c側で測定するとmyreduce2だけ30msecかかるので
    不審に思い、それまでblock=(6,8), thread=(256)で起動していたところを
    block=(1,8),thread=1024としたところ、30msecから8msecへ短縮できた。
    block=(2,8),thread=768へ変更したところ、cuda profiler上での実行時間は
    8msecから6msecになるがc側では8msecから11msecになった。
    thread生成のコストがあるようで、block=(1,8)、thread=1024が最適な模様。
    100elipsoid*100parallel*(x=256,y=512,z=153)で240msec程度になった。

ver12: 以下のメモに従い、reduceをしなくて済むように改変した。
    -> reduceにかかっていた時間が丸ごと浮いた。
    100elipsoid*100parallel*(x=256,y=512,z=153)で150msec程度になった。

    rss(1) = sum((p-m)^2)
           = sum(p^2-2*p*m+m^2)
           = sum(p^2) - 2*sum(p*m) + sum(m^2)

    m := n1+n2 として
    rss(2) = sum(p^2) - 2*sum(p*n1+p*n2) + sum(n1^2+2*n1*n2+n2^2)
           = sum(p^2) - 2*sum(p*n1) - 2*sum(p*n2) + sum(n1^2) + 2*sum(n1*n2) + sum(n2^2)
           = rss(1) - 2*sum(p*n2) + 2*sum(n1*n2) + sum(n2^2)
           = rss(1) - 2*sum((p-n1)*n2) + sum(n2^2)
           = rss(1) - sum((2*(p-n1)-n2)*n2)

    m := k1+k2+k3 とすると
    rss(3) = sum((p-((k1+k2)+k3))^2)
           = sum(p^2 - 2*p*((k1+k2)+k3) + ((k1+k2)+k3)^2)
           = sum(p^2) - 2*sum(p*((k1+k2)+k3)) + sum(((k1+k2)+k3)^2)
           = sum(p^2) - 2*sum(p*k1+p*k2) - 2*sum(p*k3) + sum((k1+k2)^2) + 2*sum((k1+k2)*k3) + sum(k3^2)
           = rss(2) - 2*sum(p*k3) + 2*sum((k1+k2)*k3) + sum(k3^2)
           = rss(2) - 2*sum((p-(k1+k2))*k3) + sum(k3^2)
           = rss(2) - sum( (2*(p-(k1+k2))-k3) * k3)

    それまでの残差2乗和と、現時点におけるその地点での残差の値が保存されていれば、k3を追加した際の残差2乗和が計算できる。
    画像サイズ分のメモリを総ナメすると8msかかるが、各elipsoidの計算には2msしかかからない。
    「現時点におけるその地点での残差の値」が保存されていれば、残差2乗和をその場で計算できる。
    メモリの初期化分のコストはかかるが、reduceのコストが丸々浮くはず。
    ただしfloatで逐次計算しているので精度は下がるはず。


ver13: MAX_NUM_ELIPSOIDが157を超える例があるようなので、elipsoidのパラメータを
    コンスタントメモリからグローバルメモリに移し、L1キャッシュを使うことにした。
    初代の頃より読み込み回数がかなり少なくなったので、所要時間はほぼ同じだった。




elipsoid_mex:
ver1: 初期型
ver2: fmath (+openmp)
ver3: fmath+xybyak (+openmp) matlabの5倍速程度に。
ver4: 2セット以上のパラメータが与えられた時の動作バグを修正
ver5: 残差計算用のthrsholdを導入
ver6: ver6: ver5の変更を戻し、cuda化を見越して並列化・単精度化
      しかしfmath::expの挙動がおかしいため標準のexpに戻した。
ver7: elipsoid_cuda_13に基づいてC版作成。double使用。並列対応あり。
    -> 100elipsoid*100parallel*(x=256,y=512,z=153)で2500msec程度。
       1parallelでも60msec程度かかる(cuda版は15msec程度)。
       並列数にかかわらず、cuda版の方が早いようだ。
       内部演算がdoubleなので、検算用に使う?
ver8: Intel compilerの自動ベクトル化+openmpを使ってAVX対応を図る。
       コードは最も簡単な形に戻し、xyz方向の範囲制約も除いた
       最大輝度値の1/2^16を超えないピクセルについては,積極的に枝刈りして計算量を減らす。
       openmpを任意のスレッド数で実行するために、im_diffなど一時記憶を撤廃し、
       ピクセル単位のループを外側に、楕円体単位のループを内側に配置しようとしたが、
       Intel compilerの自動ベクトル化は最内周のループにしか適用できず、
       ベクトル長が短いとintel VML (MKL)の速度が生かせないので断念。
       44*43*56のimgに対し楕円体3つを500個体fitさせて、平均50ms。
ver9: ver8のうち、z方向について双曲線へ拡張し、枝刈りをcutした。
ver10: 以下のメモに従って、演算回数の低減を試みた。
    31*28*55のimgに対し楕円体4つを500個体fitさせて、平均60msec(ただしsingle thread)。
    同条件でver8のとき、6threadで平均70msec -> total 420msecなので、7倍高速化したことになる。
    このとき、initで約65Mclock, loopで約60Mclock, reductionで約18Mclockだった。
    38*51*79のimgに対し楕円体4つを500個体fitさせて、平均120msec(ただしsingle thread)。
    同条件でver8は110msecだったので、約6倍高速化したことになる。
    このとき、initで約140Mclock, loopで約100Mclock, reductionで約45Mclockだった。
    init+loop部分は確実にver8より高速化したと思われる。
    reductionの所要時間はほかと比べて十分低いように見える。
ver11: mklで高速化。reduction部分をSIMD intrinsicsで記述するなど。
    38*51*79のimgに対し楕円体4つを500個体fitさせて、平均110msec(ただしsingle thread)。
    初期値のテーブルを計算するようにしてみたが、コンパイラ最適化のせいでrdtscの位置が変わったのか、
    テーブルを使わないバージョンと比べるのが困難だった。
ver12: コアループもSIMD intrinsicsで記述した。
     38*51*79のimgに対し楕円体4つを500個体fitさせて、平均100msec(ただしsingle thread)。
     38*51*79*4*500=約300M pixelなので、およそ1[clock/pixel]まで到達した。
ver13: 
    回転角が0のとき回転行列が単位行列に一致しない形式になっていたため修正。
    (過去のversionとの互換性はなくなることに注意)
    linuxのintel compilerでもcompileできるようにinclude節やdefine節を追加
ver14:
    楕円体の重複を防ぐため、罰則を導入した。
    (拡大・回転行列で楕円体中心間の距離をスケールしたものをscaled_distanceとする。
    scaled_distanceが小さいことは楕円体間の重複が大きいことを意味する。
    残差2乗和を penalty = (scaled_distance + penalty_parameter)/scaled_distance 倍することで、
    scaled_distanceがpenalty_parameterと比べて比較的小さい(楕円体が重複している)とき
    大きめのpenaltyがかかり、評価関数値が悪化する。)
ver15:
    scaled_distanceの幾何平均に基づいて罰則を計算するよう改良
ver16:
    輝度値が低いとき、scaled_distanceを計算しないよう改良
ver17:
    倍精度対応
ver18:
    残差2乗和の計算時に精度を担保するよう、加算アルゴリズムを変更(Knuth,Kahan)
ver19:
    計算した楕円体のデータを残差に直して蓄積するオプションを導入した
    楕円体を計算したボクセルでのみ残差2乗和の計算を更新するオプションを導入した
ver20:
    楕円体を計算する範囲に制限をかけられるよう変更した
ver21:
    内部で勾配を計算できるよう変更
    setofv関数内でlocal変数名を変更 scale -> log_scale
ver22:
    flag_em を 1にすると、所与のパラメータを初期値としてEm-likeな最適化を実行するよう変更


以下未定: 
svml を vmlにするため、SIMDベクター単位ではなくまとめて演算?
exp関数の精度を下げて高速化?
メモリ読み書きを考慮して、楕円体の輝度を計算した際に残差も同時に計算する?      
        


ver10におけるメモ：楕円体の輝度を求める際、指数関数は毎回計算しなくても良いのではないか？

    自然な実装：各点において、affine変換+指数関数を計算
        xtmp = xpos-xc
        xprime = xtmp*rotxx + ytmp*rotyx + ztmp*rotzx
        dist = xprime^2 + yprime^2 + zprime^2
        intensity += scale * exp(-dist)

    差分を考えた場合：
    任意の2点X1=(x,y,z)とX2=(x+dx,y+dy,z+dz)を考えた時、楕円の中心をXC=(xc,yc,zc)として、
        xpos1 = x - xc
        ypos1 = y - yc
        zpos1 = z - zc
        xpos2 = xpos1 + dx 
        ypos2 = ypos1 + dy
        zpos2 = zpos1 + dz
        
        xprime1 = xpos1*rotxx + ypos1*rotyx + zpos1*rotzx
        xprime2 = xpos2*rotxx + ypos2*rotyx + zpos2*rotzx
                = (xpos1+dx)*rotxx + (ypos1+dy)*rotyx + (zpos1+dz)*rotzx
                = xprime1 + (dx*rotxx + dy*rotyx + dz*rotzx)
        xprimed	= xrpime2 - xprime1
                = dx*rotxx + dy*rotyx + dz*rotzx
    同様に
        yprimed = dx*rotxy + dy*rotyy + dz*rotzy
        zprimed = dx*rotxz + dy*rotyz + dz*rotzz

        dist1 = xprime1^2 + yprime1^2 + zprime1^2
        dist2 = xprime2^2 + yprime2^2 + zprime2^2
              =  (xprime1 + xprimed)^2 
               + (yprime1 + yprimed)^2
               + (zprime1 + zprimed)^2
              = dist1 + 2*xprime1*xprimed + xprimed^2
                      + 2*yprime1*yprimed + yprimed^2
                      + 2*zprime1*zprimed + zprimed^2

        intensity1_add = scale * exp(-dist1)
        intensity2_add = scale * exp(-dist2)
                       = scale * exp(-dist1 - 2*xprime1*xprimed - xprimed^2
                                            - 2*yprime1*yprimed - yprimed^2
                                            - 2*zprime1*zprimed - zprimed^2)
                       = scale * exp(-dist1) 
                               * exp(- 2*xprime1*xprimed - xprimed^2
                                     - 2*yprime1*yprimed - yprimed^2
                                     - 2*zprime1*zprimed - zprimed^2)
    ここで
        dx = 0
        dy = 0
    とすると、
        xprimed = rotzx*dz
        yprimed = rotzy*dz
        zprimed = rotzz*dz
    から
        intensity2_add = scale * exp(-dist1) 
                               * exp(- 2*xprime1*(rotzx*dz) - (rotzx*dz)^2
                                     - 2*yprime1*(rotzy*dz) - (rotzy*dz)^2
                                     - 2*zprime1*(rotzz*dz) - (rotzz*dz)^2)
                       = intensity1_add
                               * exp( -2*(xprime1*rotzx+yprime1*rotzy+zprime1*rotzz)*dz
                                       - (rotzx^2+rotzy^2+rotzz^2)*dz^2 )
                       = intensity1_add * exp( - az*dz - bz*dz^2)
                       = intensity1_add * exp( - (bz*dz + az) * dz )
    となる。ただし
        az = 2*(xprime1*rotzx+yprime1*rotzy+zprime1*rotzz)
        bz = (rotzx^2+rotzy^2+rotzz^2)
    である。さらに、
        dz_(0)   = 0
        dz_(n)   = n 
        dz_(n+1) = dz_(n) + 1
    としてzを更新していくとき、
        intensity_add_(0)  = intensity1_add
                           = intensity1_add*exp(-az*dz_(0)-bz*dz_(0)^2)
                           = intensity1_add*exp(-(bz*dz_(0)+az)*dz_(0))
        intensity_add_(n)  
            = intensity1_add*exp(-(bz*dz_(n)+az)*dz_(n))
            = intensity1_add*exp(-(bz*(dz_(n-1)+1)+az)*(dz_(n-1)+1))
            = intensity1_add*exp(-(bz*dz_(n-1)+az + bz)*(dz_(n-1) + 1))
            = intensity1_add*exp(-((bz*dz_(n-1)+az)*dz_(n-1) + bz*dz_(n-1)+ (bz*dz_(n-1)+az+bz)))
            = intensity1_add*exp(-((bz*dz_(n-1)+az)*dz_(n-1) + 2*bz*dz_(n-1)+az+bz )))
            = intensity1_add*exp(-((bz*dz_(n-1)+az)*dz_(n-1)))*exp(-(2*bz*dz_(n-1)+az+bz))
            = intensity_add_(n-1)*exp(-(2*bz*dz_(n-1)+az+bz))
            = intensity_add_(n-1)*exp(-2*bz)^dz_(n-1)*exp(-(az+bz))
            = intensity_add_(n-1)*fz(n-1)
    ただし、
        fz_(n) 	= exp(-2*bz)^dz_(n)*exp(-(az+bz))
                = exp(-2*bz)^(dz_(n-1)+1)*exp(-(az+bz))
                = exp(-2*bz)^(dz_(n-1)  )*exp(-(az+bz))*exp(-2*bz)
                =             fz_(n-1)                 *exp(-2*bz)
    とした。
    
    以上の式から、各z stepにおいて、exp(-2*bz)を使って、intensity_add_(n)とfz_(n)を更新していけば、
    更新時は積算2回のみ ( = 2 [FLOPs/step])でpixelの輝度値を求めることができる。
    積算だけなので、denormalやunderflowが起きなければ誤差の蓄積もない。
    通常の演算だと、12個の積和算 + 指数関数(doubleの場合>8FLOPs, singleの場合>3FLOPs)で合計>15～20 FLOPs程度かかる。
    7～10倍高速になるかも？？
    xyのtableについて並列化すれば、AVX(256bit SIMD)でさらに4倍(double時) ～ 8倍(single時)高速になりそう。

    denormalを考えると、single floatでは不安(exp(-87)はおよそrealmin('single')と一致する)。
    z方向には輝度は単調減少か増加後減少する場合があるが、減少後増加することはないので、
    z更新中にdenormalが生じた場合は放置して(ゼロで埋めて)良い。
    したがって、z方向の計算開始点でdenormalが発生しなければ、精度は保たれる。
    
    z方向の計算開始点をどうやって求めるかが次の問題となるが、関数fz_(n)は
        fz_(n) = exp(-2*bz)^dz_(n)*exp(-(az+bz))
        dz_(n) = n
    なので、初項exp(-(az+bz))、公比exp(-2*bz)の等比数列である。第n項までの積は
        prod(fz_(k),k=0...n) = exp(-2*bz)^0*exp(-(az+bz)) * exp(-2*bz)^1*exp(-(az+bz)) * ... * exp(-2*bz)^n*exp(-(az+bz))
                             = exp(-2*bz)^(0+1+...+n)*exp(-(az+bz))^n
                             = exp(-2*bz)^(n*(n+1)/2)*exp(-(az+bz))^n
    より、
        intensity_add_(n) = intensity_add_(n-1) * fz_(n-1)
                          = intensity_add_(n-2) * fz_(n-2) * fz_(n-1)
                          = intensity_add_(0)   * fz_(0)   * fz_(1)   * ... * fz_(n-1)
                          = intensity_add_(0)   * prod(fz_(k),k=0...n-1)
                          = intensity_add_(0)   * exp(-2*bz)^((n-1)*n/2)*exp(-(az+bz))^n
                          = intensity_add_(0)   * exp(-bz)^((n-1)*n)*exp(-(az+bz))^n
                          = intensity_add_(0)   * exp(-bz)^(n^2-n)*exp(-az)^n*exp(-bz)^n
                          = intensity_add_(0)   * exp(-bz)^(n^2)*exp(-az)^n
    となる。ここで、対数をとって、    
        log(intensity_add_(n)) 
            = log(intensity_add_(0)   * exp(-bz)^(n^2)*exp(-az)^n )
            = log( exp(log(intensity_add_(0))) * exp(-bz * n^2)*exp(-az *n) )
            = log( exp(log(intensity_add_(0)) -bz*n^2 -az*n) )
            = log(intensity_add_(0)) -bz*n^2 -az*n
                
    これがdenormalしないような閾値Thr(=-87)を超えるには、
        log(intensity_add_(0)) -bz*n^2 -az*n > thr
        ->  bz*n^2 + az*n + thr - log(intensity_add_(0)) < 0
    ここでbz>0であるので、左辺はnについて下に凸な関数であるから、
        det = az^2-4*bz*(thr-log(intensity_add_(0)))
            > 0
    の時に限って実数解が存在して、解の公式より
        (-az-sqrt(det))/(2*bz) < n < (-az+sqrt(det))/(2*bz)
    といえる。

    ただし、解の公式による求解はしばしば桁落ちを起こすので注意が必要である。
    絶対値が小さい方の解は「解と係数の関係」から求めるほうが安全で、
        
        (小さい方の解) = 0次の係数 / (2次の係数 * 大きい方の解)
                       = (thr-log(intensity_add_(0))) / (az * (-az-sqrt(det))/(2*bz))

    今回はnを小数点以下切り上げ/切り捨てでfixして使うので、精度はあまり必要ではない。
    桁落ちが起きるのは、2つの解のオーダーが著しく異なる場合だが、今回はnの取りうる範囲がわかっているので、
    端点でまず条件が成立するか確認して、不成立の時だけ公式に沿って求解を実行すると良いかもしれない。
    
    上記の結論をまとめると、手順としては、
        (初期化時に実行)
        1. 計算が必要な点(x,y,z)について、z範囲を調べ、(x,y)の組をSIMDベクトルにパックする。
        2. SIMDベクトル毎にzの最大値/最小値を調べ、必要な範囲の元画像の輝度値をSIMDベクトルに格納する。
        3. あわせて、残差を計算しないvoxelについてのマスクを作成する。
        (各パラメータセットごとに実行)
        4. (x,y)の各組について、2次方程式を解いて、denormalしないようなz範囲を求める。
        5. denormalしない最小のzについて、intensity_add_(0)とfz_(0)を計算する。
        6. 更新式に従ってzを更新していく。有効なz範囲の外では0を代入するなどして対処する。



    総和計算にどのくらい時間かかっている？？
    デフォルト(/fp:precise)だと、加算順序の入替えによる最適化は結果を変化させてしまうため禁止されている(ANSI C)ので、
    スループットではなくてレイテンシが効いてきてしまう。
    /fp:fastにすると加算順序を入替えてくれるようだが、どのくらい最適化されるかは不明。
    intel C++ compilerは既定値が/fp:fastになっているようだ。
    

    ・初期化 -> 面倒なので普通に計算
    ・残差計算時のマスク -> 
    ・reduction -> 
    ・複数の楕円体の和を取ることを考えると、生成した画像を一時im_synthなどに保管するほうが良い？
    


%%%%%%%%%% 


*/


// include header file for matlab
#include "mex.h"
#include "matrix.h"

// include header file for intel mkl
//#include "mkl_vml.h"
#include "mkl.h"

// include header file for container vector
#include <vector>

// include header for OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

//// detect memory leak using msvc
//#define _crtdbg_map_alloc
//#include <stdlib.h>
//#include <crtdbg.h>

// include header file for c standard
#ifdef _MSC_VER
    #include <malloc.h>
#else
    #include <stdlib.h>
    static inline void *_aligned_malloc(size_t size, size_t alignment)
    {
        void *p;
        int ret = posix_memalign(&p, alignment, size);
        return (ret == 0) ? p : 0;
    }
    #define _aligned_free(a)    free(a)
#endif

#include <math.h>
#include <stdio.h>
#ifdef _MSC_VER
    #include <intrin.h>
#else
    #include <x86intrin.h>
#endif

/* parameter settings */
#define NUM_PARAMS_INPUT 10 // 10 pararameters per an elipsoid for input
#define NUM_PARAMS_PROC 24  // 24 pararameters per an elipsoid for calc elipsoid
#define NUM_PARAMS_CALCG 10 // 10 parameters per an elipsoid for calc_G
//#define DBL_MAX 1.7976931348623158e+308 // max limit of double float
#define M_PI 3.14159265358979323846 // pi


/* SIMD Vector settings */
//#define ENABLE_SINGLE_PRECISION // uncomment if you use single float for calculation.

#if defined(ENABLE_SINGLE_PRECISION)
    #define FLT_MAX_PRECISION 3.402823466e+38F // max limit of single float
    #define FLT_EXP_MIN_INPUT -87.0 // minimum limit of input for exponential function, single float

    #define SIMD_ELEMENT_TYPE float // single_float
    #define SIMD_VECTOR_LENGTH 8  // AVX is 256-bit SIMD (single_float*8) 
    #define __MM __m256  // AVX is 256-bit SIMD (single_float*8)
    #define __MMI __m256i // AVX is 256-bit SIMD (32bit integer*8)

    #define _MM_MASKLOAD(a,b) _mm256_maskload_ps((a), (_mm256_castps_si256((b))) )
    #define _MM_MOVEMASK(a)   _mm256_movemask_ps((a))
    #define _MM_SET1(a)       _mm256_set1_ps((a))
    #define _MM_CMP(a,b,c)    _mm256_cmp_ps((a), (b), (c))
    #define _MM_AND(a,b)      _mm256_and_ps((a), (b))
    #define _MM_OR(a,b)       _mm256_or_ps((a), (b))
    #define _MM_MAX(a,b)      _mm256_max_ps((a), (b))
    #define _MM_MIN(a,b)      _mm256_min_ps((a), (b))
    #define _MM_ADD(a,b)      _mm256_add_ps((a), (b))
    #define _MM_SUB(a,b)      _mm256_sub_ps((a), (b))
    #define _MM_MUL(a,b)      _mm256_mul_ps((a), (b))
    #define _MM_DIV(a,b)      _mm256_div_ps((a), (b))
    #define _MM_SQR(a)        _mm256_mul_ps((a), (a))
    #define _MM_SQRT(a)       _mm256_sqrt_ps((a))
    #define _MM_EXP(a)        _mm256_exp_ps((a))   // intel compiler intrinsics
    #define _MM_CEIL(a)       _mm256_ceil_ps((a)) // intel compiler intrinsics
    #define _MM_FLOOR(a)      _mm256_floor_ps((a)) // intel compiler intrinsics
    #define _MM_CVTF2I(a)     _mm256_cvtps_epi32((a))
    #define _MM_CVTI2F(a)     _mm256_cvtepi32_ps((a))
    #define _MM_BLENDV(a,b,c) _mm256_blendv_ps((a), (b), (c))
    #define _MM_PERM(a,b)     _mm256_permute_ps((a), (b))    
    #define _MM_PERM2(a,b,c)  _mm256_permute2f128_ps((a), (b), (c))
    #define _MM_TESTZ(a,b)    _mm256_testz_ps( (a), (b) )
    
    typedef union {
    uint32_T uint32;
    float f;
    } aliasing;
    aliasing aliasing_mask_true  = { 0xFFFFFFFF };
    aliasing aliasing_mask_false = { 0x00000000 };
    float mask_true  = aliasing_mask_true.f;
    float mask_false = aliasing_mask_false.f;
    __MM mask_true_mm = _MM_SET1(mask_true);
    __MM mask_false_mm = _MM_SET1(mask_false);
    #define _MM_NOT(a)   _mm256_andnot_ps((a), (mask_true_mm))
    #define _MM_TESTZ1(a) _mm256_testz_ps( (a), (mask_true_mm))
    __MM addx_mm = _mm256_setr_ps(0,1,2,3,4,5,6,7);


#else
    #define FLT_MAX_PRECISION 1.7976931348623158e+308 // max limit of double float
    #define FLT_EXP_MIN_INPUT -708.0 // minimum limit of input for exponential function, double float

    #define SIMD_ELEMENT_TYPE double // double_float
    #define SIMD_VECTOR_LENGTH 4  // AVX is 256-bit SIMD (double_float*4) 
    #define __MM __m256d  // AVX is 256-bit SIMD (double_float*4)
    #define __MMI __m128i // SSE2 is 128-bit SIMD (32bit integer*4)

    #define _MM_MASKLOAD(a,b) _mm256_maskload_pd((a), (_mm256_castpd_si256((b))) )
    #define _MM_MOVEMASK(a)   _mm256_movemask_pd((a))
    #define _MM_SET1(a)       _mm256_set1_pd((a))
    #define _MM_CMP(a,b,c)    _mm256_cmp_pd((a), (b), (c))
    #define _MM_AND(a,b)      _mm256_and_pd((a), (b))
    #define _MM_OR(a,b)       _mm256_or_pd((a), (b))
    #define _MM_MAX(a,b)      _mm256_max_pd((a), (b))
    #define _MM_MIN(a,b)      _mm256_min_pd((a), (b))
    #define _MM_ADD(a,b)      _mm256_add_pd((a), (b))
    #define _MM_SUB(a,b)      _mm256_sub_pd((a), (b))
    #define _MM_MUL(a,b)      _mm256_mul_pd((a), (b))
    #define _MM_DIV(a,b)      _mm256_div_pd((a), (b))
    #define _MM_SQR(a)        _mm256_mul_pd((a), (a))
    #define _MM_SQRT(a)       _mm256_sqrt_pd((a))
    #define _MM_EXP(a)        _mm256_exp_pd((a))   // intel compiler intrinsics
    #define _MM_CEIL(a)       _mm256_ceil_pd((a)) // intel compiler intrinsics
    #define _MM_FLOOR(a)      _mm256_floor_pd((a)) // intel compiler intrinsics
    #define _MM_CVTF2I(a)     _mm256_cvtpd_epi32((a))
    #define _MM_CVTI2F(a)     _mm256_cvtepi32_pd((a))
    #define _MM_BLENDV(a,b,c) _mm256_blendv_pd((a), (b), (c))
    #define _MM_PERM(a,b)     _mm256_permute_pd((a), (b))    
    #define _MM_PERM2(a,b,c)  _mm256_permute2f128_pd((a), (b), (c))
    #define _MM_TESTZ(a,b)    _mm256_testz_pd( (a), (b) )

    typedef union {
    uint64_T uint64;
    double f;
    } aliasing;
    aliasing aliasing_mask_true  = { 0xFFFFFFFFFFFFFFFF };
    aliasing aliasing_mask_false = { 0x0000000000000000 };
    double mask_true  = aliasing_mask_true.f;
    double mask_false = aliasing_mask_false.f;
    __MM mask_true_mm  = _MM_SET1(mask_true);
    __MM mask_false_mm = _MM_SET1(mask_false);
    #define _MM_NOT(a)   _mm256_andnot_pd((a), (mask_true_mm))
    #define _MM_TESTZ1(a) _mm256_testz_pd( (a), (mask_true_mm))
    __MM addx_mm = _mm256_setr_pd(0,1,2,3);

#endif

//#define MEMORY_ALIGNMENT 32 // alignment (bit) for memory allocation
#define MEMORY_ALIGNMENT 64 // alignment (bit) for memory allocation
typedef SIMD_ELEMENT_TYPE (*fvec)[SIMD_VECTOR_LENGTH]; // SIMD Vector
typedef               int (*ivec)[SIMD_VECTOR_LENGTH]; // SIMD Vector

#ifndef _MSC_VER
// macro of max and min for gcc
#define __max(x, y) (((x) > (y)) ? (x) : (y))
#define __min(x, y) (((x) < (y)) ? (x) : (y))
#endif


// #define MY_ASSERT 0 // use clock when defined


// declaration of static variables
SIMD_ELEMENT_TYPE *mask = NULL; // buffer for mask
SIMD_ELEMENT_TYPE *im_orig = NULL; // buffer for original image
SIMD_ELEMENT_TYPE *im_synth = NULL; // buffer for synthesized image
//SIMD_ELEMENT_TYPE *im_init_int = NULL; // table for initialization of image calculation
//SIMD_ELEMENT_TYPE *im_init_fz = NULL; // table for initialization of image calculation
//SIMD_ELEMENT_TYPE *im_end_int = NULL; // table for end point of image calculation
SIMD_ELEMENT_TYPE *idx_Zstart = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *idx_Zend = NULL; // buffer for end point of image calculation
SIMD_ELEMENT_TYPE *mem_int0 = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *mem_fz0 = NULL; // buffer for start point of image calculation
SIMD_ELEMENT_TYPE *mem_int = NULL; // buffer for image calculation
SIMD_ELEMENT_TYPE *mem_fz = NULL; // buffer for image calculation
SIMD_ELEMENT_TYPE *new_mu_x = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23 = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33 = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_x_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23_err = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33_err = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg_err = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_x_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_y_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_mu_z_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_11_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_12_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_13_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_22_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_23_tmp = NULL; // 
SIMD_ELEMENT_TYPE *new_S_33_tmp = NULL; // 
SIMD_ELEMENT_TYPE *sumrzg_tmp = NULL; // 
SIMD_ELEMENT_TYPE *params = NULL; // buffer for parameters 
SIMD_ELEMENT_TYPE *params_proc = NULL; // buffer for parameters
SIMD_ELEMENT_TYPE *volume = NULL; // volume of elipsoid, without pi
SIMD_ELEMENT_TYPE *params_calcG = NULL; // buffer for variables in calc_G

SIMD_ELEMENT_TYPE thrint, thrdist, tol;

__MM *mask_mm     = NULL;
__MM *im_orig_mm  = NULL;
__MM *im_synth_mm  = NULL;
//__MM *imInitInt_mm = NULL;
//__MM *imInitFz_mm  = NULL;
//__MM *imEndInt_mm  = NULL;
__MM *idx_Zstart_mm = NULL;
__MM *idx_Zend_mm = NULL;
__MM *mem_int0_mm = NULL;
__MM *mem_fz0_mm = NULL;
__MM *mem_int_mm = NULL;
__MM *mem_fz_mm = NULL;
__MM *new_mu_x_mm = NULL;
__MM *new_mu_y_mm = NULL;
__MM *new_mu_z_mm = NULL;
__MM *new_S_11_mm = NULL;
__MM *new_S_12_mm = NULL;
__MM *new_S_13_mm = NULL;
__MM *new_S_22_mm = NULL;
__MM *new_S_23_mm = NULL;
__MM *new_S_33_mm = NULL;
__MM *sumrzg_mm = NULL;
__MM *new_mu_x_err_mm = NULL;
__MM *new_mu_y_err_mm = NULL;
__MM *new_mu_z_err_mm = NULL;
__MM *new_S_11_err_mm = NULL;
__MM *new_S_12_err_mm = NULL;
__MM *new_S_13_err_mm = NULL;
__MM *new_S_22_err_mm = NULL;
__MM *new_S_23_err_mm = NULL;
__MM *new_S_33_err_mm = NULL;
__MM *sumrzg_err_mm = NULL;
__MM *new_mu_x_tmp_mm = NULL;
__MM *new_mu_y_tmp_mm = NULL;
__MM *new_mu_z_tmp_mm = NULL;
__MM *new_S_11_tmp_mm = NULL;
__MM *new_S_12_tmp_mm = NULL;
__MM *new_S_13_tmp_mm = NULL;
__MM *new_S_22_tmp_mm = NULL;
__MM *new_S_23_tmp_mm = NULL;
__MM *new_S_33_tmp_mm = NULL;
__MM *sumrzg_tmp_mm = NULL;
__MM *params_proc_mm = NULL;
__MM *rss_current_mm = NULL;
__MM *rss_current_err_mm = NULL;

__MM mtwos_mm = _MM_SET1(-2);
__MM mones_mm = _MM_SET1(-1);
__MM mhalf_mm = _MM_SET1(-0.5);
__MM zeros_mm = _MM_SET1(0);
__MM ones_mm  = _MM_SET1(1);
__MM twos_mm  = _MM_SET1(2);
__MM numz_mm, numzm1_mm;

double rss_orig_inv, score;

int *idx_validelip = NULL;
int *idx_minZstart = NULL;
int *idx_maxZend   = NULL;
//int *flag_invalid_allZ = NULL;

int maxLengthZ, fixmu;

unsigned int *idx_VXstart = NULL;
unsigned int *idx_Ystart = NULL;
unsigned int *idx_VXend = NULL;
unsigned int *idx_Yend = NULL;

unsigned int numx, numy, numz, padlength, numXWithPad, numVectorX, maxiter;
unsigned int numelipsoids = 0;


// The table indicates all of the voxels at a specified xy index of SIMD vector are lower than thrint.
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_mask_xy(0, std::vector<unsigned int>(0)); 

// The table contains indices of valid elipsoids at a specified xy index of SIMD vector. 
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_elipsoid_xy(0, std::vector<unsigned int>(0)); 

// The table contains indices of valid elipsoids at a specified z index of SIMD vector. 
// Will be used in mstep.
std::vector< std::vector<unsigned int> > table_elipsoid_z(0, std::vector<unsigned int>(0)); 


const mwSize *dims;
mwSize numdims;
mxArray *prhs_mldivide[2] = {NULL,NULL};
mxArray *plhs_mldivide[1] = {NULL};





/* function prototype */
void setsynth(double *psynth);
void setresult(double *presult);
void setparam_init(const mxArray *prhs0);
void setparam();
static void closefun(void);
void free_im();
void free_params();
void initialize(int nrhs, const mxArray *prhs[]);
inline void twosum(__MM *sum_mm, __MM *err_mm, __MM *add_mm);
inline void twosub(__MM *sum_mm, __MM *err_mm, __MM *sub_mm);
inline void twosum(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *addnum);
inline void twosub(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *subnum);
inline int hmax(__MM sbj_mm);
inline int hmin(__MM sbj_mm);
void em_optim();
void calc_G();
void calc_H();
void em_estep();
void em_m1step();
void em_m2step();



inline void twosum(__MM *sum_mm, __MM *err_mm, __MM *add_mm){
    __MM predictor_mm = _MM_ADD(*sum_mm,*add_mm);
    __MM corrector_mm = _MM_SUB(predictor_mm,*sum_mm);
    *err_mm = _MM_ADD(*err_mm,_MM_ADD(_MM_SUB(*sum_mm,_MM_SUB(predictor_mm,corrector_mm)),
                                      _MM_SUB(*add_mm,corrector_mm)));
    *sum_mm = predictor_mm;
}


inline void twosub(__MM *sum_mm, __MM *err_mm, __MM *sub_mm){
    __MM predictor_mm = _MM_SUB(*sum_mm,*sub_mm);
    __MM corrector_mm = _MM_SUB(predictor_mm,*sum_mm);
    *err_mm = _MM_ADD(*err_mm,_MM_SUB(_MM_SUB(*sum_mm,_MM_SUB(predictor_mm,corrector_mm)),
                                      _MM_ADD(*sub_mm,corrector_mm)));
    *sum_mm = predictor_mm;
}


inline void twosum(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *addnum){
    volatile SIMD_ELEMENT_TYPE predictor = *sumnum + *addnum;
    volatile SIMD_ELEMENT_TYPE corrector = predictor - *sumnum;
    volatile SIMD_ELEMENT_TYPE tmpnum = predictor - corrector;
    *errnum += *sumnum - tmpnum + (*addnum-corrector);
    *sumnum = predictor;
}


inline void twosub(SIMD_ELEMENT_TYPE *sumnum, SIMD_ELEMENT_TYPE *errnum, SIMD_ELEMENT_TYPE *subnum){
    volatile SIMD_ELEMENT_TYPE predictor = *sumnum - *subnum;
    volatile SIMD_ELEMENT_TYPE corrector = predictor - *sumnum;
    volatile SIMD_ELEMENT_TYPE tmpnum = predictor - corrector;
    *errnum += *sumnum - tmpnum - (*subnum+corrector);
    *sumnum = predictor;
}


inline int hmax(__MM sbj_mm) {

    __MMI sbj_mmi =_MM_CVTF2I(sbj_mm);

    int *sbj = (int *)&sbj_mmi;

    #if defined(ENABLE_SINGLE_PRECISION) // horizontal min and max for single float mm256
        return __max(__max(__max(sbj[0],sbj[1]),
                           __max(sbj[2],sbj[3])),
                     __max(__max(sbj[4],sbj[5]),
                           __max(sbj[6],sbj[7])));

    #else // horizontal min and max for double float mm256d
        return __max(__max(sbj[0],sbj[1]),
                     __max(sbj[2],sbj[3]));
    #endif

}


inline int hmin(__MM sbj_mm) {

    __MMI sbj_mmi =_MM_CVTF2I(sbj_mm);

    int *sbj = (int *)&sbj_mmi;

    #if defined(ENABLE_SINGLE_PRECISION) // horizontal min and min for single float mm256
        return __min(__min(__min(sbj[0],sbj[1]),
                           __min(sbj[2],sbj[3])),
                     __min(__min(sbj[4],sbj[5]),
                           __min(sbj[6],sbj[7])));

    #else // horizontal min and min for double float mm256d
        return __min(__min(sbj[0],sbj[1]),
                     __min(sbj[2],sbj[3]));
    #endif

}


void setsynth(double *psynth) {
    unsigned int numvxyz = numVectorX*numy*numz;


    setparam();

    for (unsigned int cvxyz=0; cvxyz<numvxyz; cvxyz++) {
        im_synth_mm[cvxyz] = zeros_mm;
    }
    
    
    // loop for elipsoid
    #ifdef _OPENMP
    unsigned int maxnumthreads = omp_get_max_threads();
    __MM *im_local_mm = (__MM *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*numvxyz*maxnumthreads,MEMORY_ALIGNMENT);
    #pragma omp parallel
    {
        unsigned int numthreads = omp_get_num_threads();
        unsigned int idx_thread = omp_get_thread_num();
        #pragma omp for
        for (unsigned int cvxyz=0; cvxyz<numvxyz*numthreads; cvxyz++) {
            im_local_mm[cvxyz] = zeros_mm;
        }

    #endif
    
    
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        unsigned int offset = celipsoid*NUM_PARAMS_PROC;

        __MM pi_k_mm    = params_proc_mm[offset +  0];
        __MM xc_mm      = params_proc_mm[offset +  1];
        __MM yc_mm      = params_proc_mm[offset +  2];
        __MM zc_mm      = params_proc_mm[offset +  3];
        __MM Sinv_11_mm = params_proc_mm[offset +  4];
        __MM Sinv_12_mm = params_proc_mm[offset +  5];
        __MM Sinv_13_mm = params_proc_mm[offset +  6];
        __MM Sinv_22_mm = params_proc_mm[offset +  7];
        __MM Sinv_23_mm = params_proc_mm[offset +  8];
        __MM Sinv_33_mm = params_proc_mm[offset +  9];
        __MM cx1_mm     = params_proc_mm[offset + 10];
        __MM cy1_mm     = params_proc_mm[offset + 11];
        __MM cx2_mm     = params_proc_mm[offset + 12];
        __MM cy2_mm     = params_proc_mm[offset + 13];
        __MM cxy_mm     = params_proc_mm[offset + 14];
        __MM cp_mm      = params_proc_mm[offset + 15];
        __MM bz_mm      = params_proc_mm[offset + 16];
        __MM exp2bz_mm  = params_proc_mm[offset + 17];

        for (unsigned int cy = idx_Ystart[celipsoid]; cy<=idx_Yend[celipsoid]; cy++) {
            __MM yd_mm  = _MM_SUB(_MM_SET1(cy), yc_mm);
            __MM yd2_mm = _MM_SQR(yd_mm);


            for (unsigned int cvx = idx_VXstart[celipsoid]; cvx<=idx_VXend[celipsoid]; cvx++) {
                __MM xd_mm = _MM_SUB(_MM_ADD(_MM_SET1(cvx*SIMD_VECTOR_LENGTH),addx_mm), xc_mm);
                __MM xd2_mm = _MM_SQR(xd_mm);
                __MM xyd_mm = _MM_MUL(xd_mm,yd_mm);

                __MM det_mm = _MM_ADD(_MM_ADD(_MM_ADD(cp_mm,_MM_MUL(cx2_mm,xd2_mm)),
                                                            _MM_MUL(cy2_mm,yd2_mm)),
                                                            _MM_MUL(cxy_mm,xyd_mm));
                __MM flagPosDet = _MM_CMP(det_mm,zeros_mm,_CMP_GE_OS);
                if (_MM_TESTZ1(flagPosDet)) { continue; }

                __MM sqrtdet_mm = _MM_SQRT(_MM_AND(det_mm,flagPosDet));
                __MM cx1xcy1y_mm = _MM_ADD(_MM_ADD(_MM_MUL(cx1_mm,xd_mm),_MM_MUL(cy1_mm,yd_mm)),zc_mm);
                __MM tmpIdxMinZ_mm =  _MM_CEIL(_MM_SUB(cx1xcy1y_mm,sqrtdet_mm)); 
                __MM tmpIdxMaxZ_mm = _MM_FLOOR(_MM_ADD(cx1xcy1y_mm,sqrtdet_mm)); 
            
                __MM zstart_mm = _MM_MAX(tmpIdxMinZ_mm, zeros_mm); 
                __MM   zend_mm = _MM_MIN(tmpIdxMaxZ_mm,numzm1_mm); 

                __MM flagInRange = _MM_AND(flagPosDet,_MM_AND(
                        _MM_CMP(zstart_mm,numzm1_mm,_CMP_LE_OS), 
                        _MM_CMP(  zend_mm, zeros_mm,_CMP_GE_OS)));
                zstart_mm = _MM_BLENDV( numz_mm,zstart_mm,flagInRange);
                //zend_mm   = _MM_BLENDV(mones_mm,  zend_mm,flagInRange);
                zend_mm = _MM_BLENDV(zeros_mm,_MM_ADD(zend_mm,ones_mm),flagInRange);

                __MM zd_mm = _MM_SUB(zstart_mm,zc_mm);

                __MM xs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                                             _MM_MUL(yd_mm,Sinv_12_mm)),
                                             _MM_MUL(zd_mm,Sinv_13_mm));
                __MM ys_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                                             _MM_MUL(yd_mm,Sinv_22_mm)),
                                             _MM_MUL(zd_mm,Sinv_23_mm));
                __MM zs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                                             _MM_MUL(yd_mm,Sinv_23_mm)),
                                             _MM_MUL(zd_mm,Sinv_33_mm));
            
                __MM intensity0_mm = _MM_EXP(_MM_MUL(mhalf_mm,_MM_ADD(_MM_ADD(_MM_MUL(xd_mm,xs_mm),
                                                                              _MM_MUL(yd_mm,ys_mm)),
                                                                              _MM_MUL(zd_mm,zs_mm))));
                __MM az_mm = zs_mm;
                __MM fz0_mm = _MM_EXP(_MM_SUB(bz_mm,az_mm));


                //int minZStart, maxZEnd;
                //hmin(&minZStart,zstart_mm);
                //hmax(  &maxZEnd,  zend_mm);
                int minZStart = hmin(zstart_mm);
                int   maxZEnd = hmax(  zend_mm);

                ////__MMI zstart_mmi =_MM_CVTF2I(zstart_mm);
                ////__MMI   zend_mmi =_MM_CVTF2I(  zend_mm);
                ////int *zstart = (int *)&zstart_mmi;
                ////int   *zend = (int *)&zend_mmi;

                ////#if defined(ENABLE_SINGLE_PRECISION) // horizontal min and max for single float mm256
                ////    int minZStart = __min(__min(__min(zstart[0],zstart[1]),
                ////                                __min(zstart[2],zstart[3])),
                ////                          __min(__min(zstart[4],zstart[5]),
                ////                                __min(zstart[6],zstart[7])));
                ////    int maxZEnd = __max(__max(__max(zend[0],zend[1]),
                ////                              __max(zend[2],zend[3])),
                ////                        __max(__max(zend[4],zend[5]),
                ////                              __max(zend[6],zend[7])));
                ////#else // horizontal min and max for double float mm256d
                ////    int minZStart = __min(__min(zstart[0],zstart[1]),
                ////                          __min(zstart[2],zstart[3]));                                            
                ////    int maxZEnd = __max(__max(zend[0],zend[1]),
                ////                        __max(zend[2],zend[3]));
                ////#endif


                //// /* set the initial values to the tables */
                ////SIMD_ELEMENT_TYPE *intensity0 = (SIMD_ELEMENT_TYPE *)&intensity0_mm;
                ////SIMD_ELEMENT_TYPE *fz0        = (SIMD_ELEMENT_TYPE *)&fz0_mm;
                ////for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
                ////    im_init_int[ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = intensity0[counter];
                ////    im_init_fz [ zstart[counter] *SIMD_VECTOR_LENGTH+counter] =        fz0[counter];
                ////    im_end_int [(1+zend[counter])*SIMD_VECTOR_LENGTH+counter] = mask_false;                    
                ////}

                /* loop for z */
                __MM intensity_mm = zeros_mm;
                __MM fz_mm = zeros_mm;
                unsigned int cvxyz;
                //for (int cz = minZStart; cz<=maxZEnd; cz++ ){
                for (int cz = minZStart; cz<maxZEnd; cz++ ){
                    /* z方向の更新式 */
                    __MM posZ_mm = _MM_SET1(cz); 
                    __MM flagstart_mm = _MM_CMP(zstart_mm,posZ_mm,_CMP_EQ_OS);
                    __MM flagend_mm   = _MM_CMP(  zend_mm,posZ_mm,_CMP_NEQ_OS);
                
                    intensity_mm = _MM_AND(_MM_ADD(_MM_MUL(intensity_mm,fz_mm),_MM_AND(intensity0_mm,flagstart_mm)),flagend_mm);
                    fz_mm        =         _MM_ADD(_MM_MUL(fz_mm,exp2bz_mm),   _MM_AND(fz0_mm,flagstart_mm));


                    //intensity_mm = _MM_AND(_MM_ADD(_MM_MUL(intensity_mm,fz_mm),imInitInt_mm[cz]),imEndInt_mm[cz]);
                    //fz_mm        =         _MM_ADD(_MM_MUL(fz_mm,exp2bz_mm),   imInitFz_mm[cz]);
                
                    cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                    #ifdef _OPENMP
                    unsigned int offset = numvxyz*idx_thread;
                    im_local_mm[offset+cvxyz] = _MM_AND(_MM_ADD(im_local_mm[offset+cvxyz],_MM_MUL(pi_k_mm,intensity_mm)),mask_mm[cvxyz]);
                    #else
                    im_synth_mm[cvxyz] = _MM_AND(_MM_ADD(im_synth_mm[cvxyz],_MM_MUL(pi_k_mm,intensity_mm)),mask_mm[cvxyz]);
                    #endif
                }

                /////* reset initial value tables */
                ////for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
                ////    im_init_int[ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = 0.0;
                ////    im_init_fz [ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = 0.0;
                ////    im_end_int [(1+zend[counter])*SIMD_VECTOR_LENGTH+counter] = mask_true;                    
                ////}

            }
        }
    } // end of omp_for

    //#pragma omp barrier

    #ifdef _OPENMP
    // aggregate the local copies into global array
    #pragma omp for
    for (unsigned int cvxyz=0; cvxyz<numvxyz; cvxyz++) {
        for (unsigned int p=0; p<numthreads; p++) {
            //im_synth_mm[cvxyz] = _MM_ADD(im_synth_mm[cvxyz],plocal[p][cvxyz]);
            unsigned int offset = numvxyz*p;
            im_synth_mm[cvxyz] = _MM_ADD(im_synth_mm[cvxyz],im_local_mm[offset+cvxyz]);
        }
    }

    ////_aligned_free(im_local_mm);
    ////im_local_mm = NULL;

    } // end of parallel

    _aligned_free(im_local_mm);
    im_local_mm = NULL;
    #endif



    /* set data */
    unsigned int cvxyz, cout;

    for (unsigned int cy=0; cy<numy; cy++) {        
        
        // for elements without padding        
        for (unsigned int cvx=0; cvx<numVectorX-1; cvx++) {
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cout  = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                    for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH; ce++) {
                    psynth[ cout + ce ] = (double) im_synth[ cvxyz*SIMD_VECTOR_LENGTH + ce ];
                }
            }
        }

        // for elements with padding
        unsigned int cvx = numVectorX - 1; {
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cout   = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
                    psynth[ cout + ce ] = (double) im_synth[ cvxyz*SIMD_VECTOR_LENGTH + ce ];                  
                }
            }
        }
    }


}


void setresult(double *presult) {

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
        presult[NUM_PARAMS_INPUT*celipsoid + 0] = params[NUM_PARAMS_INPUT*celipsoid + 0]; // pi
        presult[NUM_PARAMS_INPUT*celipsoid + 1] = params[NUM_PARAMS_INPUT*celipsoid + 1] + 1.0; // xc
        presult[NUM_PARAMS_INPUT*celipsoid + 2] = params[NUM_PARAMS_INPUT*celipsoid + 2] + 1.0; // yc
        presult[NUM_PARAMS_INPUT*celipsoid + 3] = params[NUM_PARAMS_INPUT*celipsoid + 3] + 1.0; // zc
        presult[NUM_PARAMS_INPUT*celipsoid + 4] = params[NUM_PARAMS_INPUT*celipsoid + 4]; // S_11
        presult[NUM_PARAMS_INPUT*celipsoid + 5] = params[NUM_PARAMS_INPUT*celipsoid + 5]; // S_12
        presult[NUM_PARAMS_INPUT*celipsoid + 6] = params[NUM_PARAMS_INPUT*celipsoid + 6]; // S_13
        presult[NUM_PARAMS_INPUT*celipsoid + 7] = params[NUM_PARAMS_INPUT*celipsoid + 7]; // S_22
        presult[NUM_PARAMS_INPUT*celipsoid + 8] = params[NUM_PARAMS_INPUT*celipsoid + 8]; // S_23
        presult[NUM_PARAMS_INPUT*celipsoid + 9] = params[NUM_PARAMS_INPUT*celipsoid + 9]; // S_33
    }

}


void em_optim() {

    #ifdef MY_ASSERT
        unsigned __int64 t_start=0,t_setparam=0,t_estep=0,t_m1step=0,t_m2step=0;
    #endif

    double oldscore = DBL_MAX;

    for (unsigned int citer=0; citer<maxiter; citer++) { 

        mexPrintf("citer = %d,\t", citer);
        mexEvalString("drawnow");

        if (!fixmu) {

            #ifdef MY_ASSERT
                t_start = __rdtsc();
            #endif

            setparam();
            
            #ifdef MY_ASSERT
                t_setparam += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif

            em_estep();

            #ifdef MY_ASSERT
                t_estep += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif

            em_m1step();

            #ifdef MY_ASSERT
                t_m1step += __rdtsc() - t_start;
            #endif
        }

        #ifdef MY_ASSERT
            t_start = __rdtsc();
        #endif
        
        setparam();

        #ifdef MY_ASSERT
            t_setparam += __rdtsc() - t_start;
            t_start = __rdtsc();
        #endif
        
        em_estep();

        #ifdef MY_ASSERT
            t_estep += __rdtsc() - t_start;
            t_start = __rdtsc();
        #endif

        em_m2step();

        #ifdef MY_ASSERT
            t_m2step += __rdtsc() - t_start;
        #endif

        mexPrintf("score = %15.12f\n", score);

        if ((oldscore-score)<tol) {break;}
        oldscore = score;
        

    }
    
    #ifdef MY_ASSERT
        mexPrintf("em_optim.setparam: %I64d[clocks]\n", t_setparam);
        mexPrintf("em_optim.estep: %I64d[clocks]\n", t_estep);
        mexPrintf("em_optim.m1step: %I64d[clocks]\n", t_m1step);
        mexPrintf("em_optim.m2step: %I64d[clocks]\n", t_m2step);
    #endif

}


void em_estep() {
/* update pi_k (denoted as p)*/
    
    #ifdef MY_ASSERT    
        unsigned __int64 t_start, t_memory=0, t_calcG=0, t_calcH=0, t_mldivide=0;   
    #endif

    // calc new pi
    #ifdef MY_ASSERT
        t_start = __rdtsc();
    #endif
    
    prhs_mldivide[0] = mxCreateDoubleMatrix(numelipsoids,numelipsoids,mxREAL);
    prhs_mldivide[1] = mxCreateDoubleMatrix(numelipsoids,           1,mxREAL);
    
    #ifdef MY_ASSERT
        t_memory += __rdtsc() - t_start;
        t_start = __rdtsc();
    #endif

    calc_G();

    #ifdef MY_ASSERT
        t_calcG += __rdtsc() - t_start;
        t_start = __rdtsc();
    #endif

    calc_H();

    #ifdef MY_ASSERT
        t_calcH += __rdtsc() - t_start;
        t_start = __rdtsc();
    #endif

    mexCallMATLAB(1,plhs_mldivide,2,prhs_mldivide,"mldivide");

    #ifdef MY_ASSERT
        t_mldivide += __rdtsc() - t_start;
    #endif
	
    // update pi
    double *pin = mxGetPr(plhs_mldivide[0]);
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
        params[NUM_PARAMS_INPUT*celipsoid] = (SIMD_ELEMENT_TYPE) pin[celipsoid]; // pi
        params_proc_mm[NUM_PARAMS_PROC*celipsoid] = _MM_SET1(pin[celipsoid]);
    }

    #ifdef MY_ASSERT
        t_start = __rdtsc();
    #endif

    mxDestroyArray(prhs_mldivide[0]);
    mxDestroyArray(prhs_mldivide[1]);
    mxDestroyArray(plhs_mldivide[0]);
    prhs_mldivide[0] = NULL;
    prhs_mldivide[1] = NULL;
    plhs_mldivide[0] = NULL;

    #ifdef MY_ASSERT
        t_memory += __rdtsc() - t_start;
    #endif
    
    #ifdef MY_ASSERT
        mexPrintf("em_optim.estep.t_memory: %I64d[clocks]\n", t_memory);
        mexPrintf("em_optim.estep.t_calcG: %I64d[clocks]\n", t_calcG);
        mexPrintf("em_optim.estep.t_calcH: %I64d[clocks]\n", t_calcH);
        mexPrintf("em_optim.estep.t_mldivide: %I64d[clocks]\n", t_mldivide);
    #endif

}


void em_m1step() {
/* update mean vector (denoted as mu_k) */
    
    #ifdef MY_ASSERT    
        unsigned __int64 t_start, t_init=0, t_loop=0, t_update=0;
    #endif
   
    unsigned int numthreads = 1;
    unsigned int idx_thread = 0;
    #ifdef _OPENMP
        #pragma omp parallel firstprivate(table_elipsoid_z,idx_thread)
        { 
            numthreads = omp_get_num_threads();
            idx_thread = omp_get_thread_num();
            rss_current_mm[idx_thread] = zeros_mm;
            rss_current_err_mm[idx_thread] = zeros_mm;
    #endif

    unsigned int omp_offset = numelipsoids*idx_thread;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int celipsoid=0; celipsoid<numelipsoids*numthreads; celipsoid++) {
        new_mu_x_mm[celipsoid] = zeros_mm;
        new_mu_y_mm[celipsoid] = zeros_mm;
        new_mu_z_mm[celipsoid] = zeros_mm;
        sumrzg_mm[celipsoid]   = zeros_mm;
        new_mu_x_err_mm[celipsoid] = zeros_mm;
        new_mu_y_err_mm[celipsoid] = zeros_mm;
        new_mu_z_err_mm[celipsoid] = zeros_mm;
        sumrzg_err_mm[celipsoid]   = zeros_mm;
    } // end of omp_for 

    // loop for xy vectors
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int cy = 0; cy<numy; cy++) {
        __MM posY_mm = _MM_SET1(cy);

        //for (unsigned int cvx = 0; cvx<numVectorX; cvx++) {
        for (unsigned int cx=0; cx<table_mask_xy[cy].size(); cx++) {
            unsigned int cvx = table_mask_xy[cy][cx];

            //if (flag_invalid_allZ[cy*numVectorX+cvx]) { continue; }

            __MM posX_mm = _MM_ADD(_MM_SET1(cvx*SIMD_VECTOR_LENGTH),addx_mm);

            unsigned int cvalid = omp_offset;
            unsigned int tableidx = cy*numVectorX+cvx;
            __MM min_zstart_mm = numz_mm;
            __MM max_zend_mm = zeros_mm;

            for (unsigned int cz=0; cz<table_elipsoid_z.size(); cz++) {
                table_elipsoid_z[cz].clear();
            }

            #ifdef MY_ASSERT
                t_start = __rdtsc();
            #endif

            // loop for elipsoid
            for (unsigned int c=0; c<table_elipsoid_xy[tableidx].size(); c++) {

                unsigned int offset = table_elipsoid_xy[tableidx][c]*NUM_PARAMS_PROC;

                __MM xc_mm      = params_proc_mm[offset +  1];
                __MM yc_mm      = params_proc_mm[offset +  2];
                __MM zc_mm      = params_proc_mm[offset +  3];
                __MM Sinv_11_mm = params_proc_mm[offset +  4];
                __MM Sinv_12_mm = params_proc_mm[offset +  5];
                __MM Sinv_13_mm = params_proc_mm[offset +  6];
                __MM Sinv_22_mm = params_proc_mm[offset +  7];
                __MM Sinv_23_mm = params_proc_mm[offset +  8];
                __MM Sinv_33_mm = params_proc_mm[offset +  9];
                __MM cx1_mm     = params_proc_mm[offset + 10];
                __MM cy1_mm     = params_proc_mm[offset + 11];
                __MM cx2_mm     = params_proc_mm[offset + 12];
                __MM cy2_mm     = params_proc_mm[offset + 13];
                __MM cxy_mm     = params_proc_mm[offset + 14];
                __MM cp_mm      = params_proc_mm[offset + 15];
                __MM bz_mm      = params_proc_mm[offset + 16];



                __MM xd_mm  = _MM_SUB(posX_mm, xc_mm);
                __MM yd_mm  = _MM_SUB(posY_mm, yc_mm);
                __MM xd2_mm = _MM_SQR(xd_mm);
                __MM yd2_mm = _MM_SQR(yd_mm);
                __MM xyd_mm = _MM_MUL(xd_mm,yd_mm);
        
                __MM det_mm = _MM_ADD(_MM_ADD(_MM_ADD(cp_mm,_MM_MUL(cx2_mm,xd2_mm)),
                                                            _MM_MUL(cy2_mm,yd2_mm)),
                                                            _MM_MUL(cxy_mm,xyd_mm));
                __MM flagPosDet = _MM_CMP(det_mm,zeros_mm,_CMP_GE_OS);
                if (_MM_TESTZ1(flagPosDet)) { continue; }

                __MM sqrtdet_mm = _MM_SQRT(_MM_AND(det_mm,flagPosDet));
                __MM cx1xcy1y_mm = _MM_ADD(_MM_ADD(_MM_MUL(cx1_mm,xd_mm),_MM_MUL(cy1_mm,yd_mm)),zc_mm);
                __MM tmpIdxMinZ_mm =  _MM_CEIL(_MM_SUB(cx1xcy1y_mm,sqrtdet_mm)); 
                __MM tmpIdxMaxZ_mm = _MM_FLOOR(_MM_ADD(cx1xcy1y_mm,sqrtdet_mm)); 
            
                __MM zstart_mm = _MM_MAX(tmpIdxMinZ_mm, zeros_mm); 
                __MM   zend_mm = _MM_MIN(tmpIdxMaxZ_mm,numzm1_mm); 

                __MM flagInRange = _MM_AND(flagPosDet,_MM_AND(
                        _MM_CMP(zstart_mm,numzm1_mm,_CMP_LE_OS), 
                        _MM_CMP(  zend_mm, zeros_mm,_CMP_GE_OS)));
                zstart_mm = _MM_BLENDV( numz_mm,zstart_mm,flagInRange);
                //zend_mm   = _MM_BLENDV(mones_mm,  zend_mm,flagInRange);
                zend_mm = _MM_BLENDV(zeros_mm,_MM_ADD(zend_mm,ones_mm),flagInRange);

                __MM zd_mm = _MM_SUB(zstart_mm,zc_mm);

                __MM xs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                                             _MM_MUL(yd_mm,Sinv_12_mm)),
                                             _MM_MUL(zd_mm,Sinv_13_mm));
                __MM ys_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                                             _MM_MUL(yd_mm,Sinv_22_mm)),
                                             _MM_MUL(zd_mm,Sinv_23_mm));
                __MM zs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                                             _MM_MUL(yd_mm,Sinv_23_mm)),
                                             _MM_MUL(zd_mm,Sinv_33_mm));
            
                __MM intensity0_mm = _MM_EXP(_MM_MUL(mhalf_mm,_MM_ADD(_MM_ADD(_MM_MUL(xd_mm,xs_mm),
                                                                              _MM_MUL(yd_mm,ys_mm)),
                                                                              _MM_MUL(zd_mm,zs_mm))));
                __MM az_mm = zs_mm;
                __MM fz0_mm = _MM_EXP(_MM_SUB(bz_mm,az_mm));


                //int tmpzmin,tmpzmax;
                //hmin(&tmpzmin,zstart_mm);
                //hmax(&tmpzmax,  zend_mm);
                int tmpzmin = hmin(zstart_mm);
                int tmpzmax = hmax(  zend_mm);
                for (int cz=tmpzmin; cz<tmpzmax; cz++) {
                    table_elipsoid_z[cz].push_back(cvalid);
                }               
                
                min_zstart_mm = _MM_MIN(min_zstart_mm,zstart_mm);
                max_zend_mm = _MM_MAX( max_zend_mm, zend_mm);
                idx_validelip[cvalid] = table_elipsoid_xy[tableidx][c];
                mem_int0_mm[cvalid] = intensity0_mm;
                mem_fz0_mm[cvalid] = fz0_mm;
                mem_int_mm[cvalid] = zeros_mm;
                mem_fz_mm[cvalid] = zeros_mm;
                idx_Zstart_mm[cvalid] = zstart_mm;
                idx_Zend_mm[cvalid] = zend_mm;
                cvalid++;           

            }
            unsigned int numvalid = cvalid;


            //int minZStart, maxZEnd;
            //hmin(&minZStart,min_zstart_mm);
            //hmax(  &maxZEnd,  max_zend_mm);
            int minZStart = hmin(min_zstart_mm);
            int   maxZEnd = hmax(  max_zend_mm);


            #ifdef MY_ASSERT
                t_init += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif


             /* loop for z */
            unsigned int cvxyz;
            for (int cz = minZStart; cz<maxZEnd; cz++ ){
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;

                if (table_elipsoid_z[cz].size()==0) { continue; }

                // update intensity of each elipsoid at the voxel
                __MM posZ_mm = _MM_SET1(cz);
                //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM flagstart_mm = _MM_CMP(idx_Zstart_mm[cvalid],posZ_mm,_CMP_EQ_OS);
                    __MM flagend_mm   = _MM_CMP(  idx_Zend_mm[cvalid],posZ_mm,_CMP_NEQ_OS);
                    __MM exp2bz_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 17];
                
                    mem_int_mm[cvalid] = _MM_AND(_MM_ADD(_MM_MUL(mem_int_mm[cvalid],mem_fz_mm[cvalid]),
                                                         _MM_AND(mem_int0_mm[cvalid],flagstart_mm)),flagend_mm);
                    mem_fz_mm[cvalid]  =         _MM_ADD(_MM_MUL(mem_fz_mm[cvalid],exp2bz_mm),
                                                         _MM_AND(mem_fz0_mm[cvalid],flagstart_mm));
                }

                // when the voxel was masked, skip the following calculations
                if (_MM_TESTZ1(mask_mm[cvxyz])) { continue; }
                
                // calc integrals 
                __MM intensity_mm = zeros_mm;
                //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM pi_k_mm   = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  0];
                    __MM tmpint_mm = _MM_MUL(mem_int_mm[cvalid],pi_k_mm);
                    intensity_mm = _MM_ADD(intensity_mm,tmpint_mm);
                
                    __MM tmpint2_mm = _MM_MUL(mem_int_mm[cvalid],tmpint_mm);
                    new_mu_x_tmp_mm[cvalid] = _MM_MUL(posX_mm,tmpint2_mm);
                    new_mu_y_tmp_mm[cvalid] = _MM_MUL(posY_mm,tmpint2_mm);
                    new_mu_z_tmp_mm[cvalid] = _MM_MUL(posZ_mm,tmpint2_mm);
                    sumrzg_tmp_mm  [cvalid] =                 tmpint2_mm;

                }
                
                __MM residual_mm = _MM_SUB(im_orig_mm[cvxyz], intensity_mm);
                __MM rss_old_mm  = _MM_SQR(im_orig_mm[cvxyz]);
                __MM rss_tmp_mm  = _MM_SQR(residual_mm);
                twosub(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_old_mm);
                twosum(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_tmp_mm);

                __MM tmpratio_mm = _MM_DIV(residual_mm,intensity_mm);
                tmpratio_mm = _MM_AND(tmpratio_mm,_MM_CMP(intensity_mm,zeros_mm,_CMP_NEQ_OQ));
            
                //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    unsigned int tmpidx = idx_validelip[cvalid] + omp_offset;
                    __MM tmp_mu_x_tmp_mm   = _MM_MUL(new_mu_x_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_mu_y_tmp_mm   = _MM_MUL(new_mu_y_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_mu_z_tmp_mm   = _MM_MUL(new_mu_z_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_sumrzg_tmp_mm = _MM_MUL(  sumrzg_tmp_mm[cvalid],tmpratio_mm);
                    twosum(&new_mu_x_mm[tmpidx],&new_mu_x_err_mm[tmpidx],  &tmp_mu_x_tmp_mm);
                    twosum(&new_mu_y_mm[tmpidx],&new_mu_y_err_mm[tmpidx],  &tmp_mu_y_tmp_mm);
                    twosum(&new_mu_z_mm[tmpidx],&new_mu_z_err_mm[tmpidx],  &tmp_mu_z_tmp_mm);
                    twosum(  &sumrzg_mm[tmpidx],  &sumrzg_err_mm[tmpidx],&tmp_sumrzg_tmp_mm);
                }
            }

            #ifdef MY_ASSERT
                t_loop += __rdtsc() - t_start;
            #endif

        }
    } // end of omp_for

    #ifdef _OPENMP

    // aggregate the local copies into global array
    #pragma omp for
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
        for (unsigned int p=1; p<numthreads; p++) {
            twosum(&new_mu_x_mm[celipsoid],&new_mu_x_err_mm[celipsoid],&new_mu_x_mm[numelipsoids*p+celipsoid]);
            twosum(&new_mu_y_mm[celipsoid],&new_mu_y_err_mm[celipsoid],&new_mu_y_mm[numelipsoids*p+celipsoid]);
            twosum(&new_mu_z_mm[celipsoid],&new_mu_z_err_mm[celipsoid],&new_mu_z_mm[numelipsoids*p+celipsoid]);
            twosum(  &sumrzg_mm[celipsoid],  &sumrzg_err_mm[celipsoid],  &sumrzg_mm[numelipsoids*p+celipsoid]);
            twosum(&new_mu_x_mm[celipsoid],&new_mu_x_err_mm[celipsoid],&new_mu_x_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_mu_y_mm[celipsoid],&new_mu_y_err_mm[celipsoid],&new_mu_y_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_mu_z_mm[celipsoid],&new_mu_z_err_mm[celipsoid],&new_mu_z_err_mm[numelipsoids*p+celipsoid]);
            twosum(  &sumrzg_mm[celipsoid],  &sumrzg_err_mm[celipsoid],  &sumrzg_err_mm[numelipsoids*p+celipsoid]);
        }
    } // end of omp_for


    } // end of parallel

    // aggregate the local copies of rss_current_mm into global array
    for (unsigned int p=1; p<numthreads; p++) {
        twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_mm[p]);
        twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_err_mm[p]);
    }
    #endif
       

    #ifdef MY_ASSERT
        t_start = __rdtsc();
    #endif

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        /* sum up new_mu_* and set to tmpnum */
        SIMD_ELEMENT_TYPE tmpnum_x     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_y     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_z     = 0.0;
        SIMD_ELEMENT_TYPE tmpden       = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_x_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_y_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_z_err = 0.0;
        SIMD_ELEMENT_TYPE tmpden_err   = 0.0;
        for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
            twosum(&tmpnum_x,&tmpnum_x_err,&new_mu_x[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_y,&tmpnum_y_err,&new_mu_y[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_z,&tmpnum_z_err,&new_mu_z[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpden,  &tmpden_err,  &sumrzg[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_x,&tmpnum_x_err,&new_mu_x_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_y,&tmpnum_y_err,&new_mu_y_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_z,&tmpnum_z_err,&new_mu_z_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpden,  &tmpden_err,  &sumrzg_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
        }
        tmpnum_x += tmpnum_x_err;
        tmpnum_y += tmpnum_y_err;
        tmpnum_z += tmpnum_z_err;
        tmpden   += tmpden_err;

        /* update mu */
        SIMD_ELEMENT_TYPE pi_k     = params[NUM_PARAMS_INPUT*celipsoid + 0];
        SIMD_ELEMENT_TYPE old_mu_x = params[NUM_PARAMS_INPUT*celipsoid + 1];
        SIMD_ELEMENT_TYPE old_mu_y = params[NUM_PARAMS_INPUT*celipsoid + 2];
        SIMD_ELEMENT_TYPE old_mu_z = params[NUM_PARAMS_INPUT*celipsoid + 3];
        
        if (pi_k==0) { continue; }

        tmpnum_x += old_mu_x*pi_k*volume[celipsoid];
        tmpnum_y += old_mu_y*pi_k*volume[celipsoid];
        tmpnum_z += old_mu_z*pi_k*volume[celipsoid];
        tmpden   +=          pi_k*volume[celipsoid];

        SIMD_ELEMENT_TYPE tmpden_inv = 1/tmpden;
        
        params[NUM_PARAMS_INPUT*celipsoid + 1] = tmpnum_x * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 2] = tmpnum_y * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 3] = tmpnum_z * tmpden_inv;

    }


    // update score
    SIMD_ELEMENT_TYPE tmpsum     = 0.0;
    SIMD_ELEMENT_TYPE tmpsum_err = 0.0;
    SIMD_ELEMENT_TYPE *rss_current     = (SIMD_ELEMENT_TYPE *)&rss_current_mm[0];
    SIMD_ELEMENT_TYPE *rss_current_err = (SIMD_ELEMENT_TYPE *)&rss_current_err_mm[0];
    for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
        twosum(&tmpsum,&tmpsum_err,&rss_current[counter]);
        twosum(&tmpsum,&tmpsum_err,&rss_current_err[counter]);
    }
    tmpsum += tmpsum_err;
    score = 1 + ((double)tmpsum) * rss_orig_inv;

    #ifdef MY_ASSERT
        t_update = __rdtsc() - t_start;
    #endif

    #ifdef MY_ASSERT    
        mexPrintf("em_optim.m1step.t_init: %I64d[clocks]\n", t_init);
        mexPrintf("em_optim.m1step.t_loop: %I64d[clocks]\n", t_loop);
        mexPrintf("em_optim.m1step.t_update: %I64d[clocks]\n", t_update);
    #endif

}


void em_m2step() {
/* update sigma vector (denoted as cv_k) */

    #ifdef MY_ASSERT    
        unsigned __int64 t_start, t_init=0, t_loop=0, t_update=0;
    #endif

    unsigned int numthreads = 1;
    unsigned int idx_thread = 0;
    #ifdef _OPENMP
        #pragma omp parallel firstprivate(table_elipsoid_z,idx_thread)
        { 
            numthreads = omp_get_num_threads();
            idx_thread = omp_get_thread_num();
            rss_current_mm[idx_thread] = zeros_mm;
            rss_current_err_mm[idx_thread] = zeros_mm;
    #endif

    unsigned int omp_offset = numelipsoids*idx_thread;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int celipsoid=0; celipsoid<numelipsoids*numthreads; celipsoid++) {
        new_S_11_mm[celipsoid] = zeros_mm;
        new_S_12_mm[celipsoid] = zeros_mm;
        new_S_13_mm[celipsoid] = zeros_mm;
        new_S_22_mm[celipsoid] = zeros_mm;
        new_S_23_mm[celipsoid] = zeros_mm;
        new_S_33_mm[celipsoid] = zeros_mm;
        sumrzg_mm[celipsoid]   = zeros_mm;
        new_S_11_err_mm[celipsoid] = zeros_mm;
        new_S_12_err_mm[celipsoid] = zeros_mm;
        new_S_13_err_mm[celipsoid] = zeros_mm;
        new_S_22_err_mm[celipsoid] = zeros_mm;
        new_S_23_err_mm[celipsoid] = zeros_mm;
        new_S_33_err_mm[celipsoid] = zeros_mm;
        sumrzg_err_mm[celipsoid]   = zeros_mm;
    } // end of omp_for 

    // loop for xy vectors
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for (unsigned int cy = 0; cy<numy; cy++) {
        __MM posY_mm = _MM_SET1(cy);

        //for (unsigned int cvx = 0; cvx<numVectorX; cvx++) {
        for (unsigned int cx=0; cx<table_mask_xy[cy].size(); cx++) {
            unsigned int cvx = table_mask_xy[cy][cx];

            //if (flag_invalid_allZ[cy*numVectorX+cvx]) { continue; }

            __MM posX_mm = _MM_ADD(_MM_SET1(cvx*SIMD_VECTOR_LENGTH),addx_mm);

            unsigned int cvalid = omp_offset;
            unsigned int tableidx = cy*numVectorX+cvx;
            __MM min_zstart_mm = numz_mm;
            __MM max_zend_mm = zeros_mm;

            for (unsigned int cz=0; cz<table_elipsoid_z.size(); cz++) {
                table_elipsoid_z[cz].clear();
            }

            #ifdef MY_ASSERT
                t_start = __rdtsc();
            #endif

            // loop for elipsoid
            for (unsigned int c=0; c<table_elipsoid_xy[tableidx].size(); c++) {

                unsigned int offset = table_elipsoid_xy[tableidx][c]*NUM_PARAMS_PROC;

                __MM xc_mm      = params_proc_mm[offset +  1];
                __MM yc_mm      = params_proc_mm[offset +  2];
                __MM zc_mm      = params_proc_mm[offset +  3];
                __MM Sinv_11_mm = params_proc_mm[offset +  4];
                __MM Sinv_12_mm = params_proc_mm[offset +  5];
                __MM Sinv_13_mm = params_proc_mm[offset +  6];
                __MM Sinv_22_mm = params_proc_mm[offset +  7];
                __MM Sinv_23_mm = params_proc_mm[offset +  8];
                __MM Sinv_33_mm = params_proc_mm[offset +  9];
                __MM cx1_mm     = params_proc_mm[offset + 10];
                __MM cy1_mm     = params_proc_mm[offset + 11];
                __MM cx2_mm     = params_proc_mm[offset + 12];
                __MM cy2_mm     = params_proc_mm[offset + 13];
                __MM cxy_mm     = params_proc_mm[offset + 14];
                __MM cp_mm      = params_proc_mm[offset + 15];
                __MM bz_mm      = params_proc_mm[offset + 16];



                __MM xd_mm  = _MM_SUB(posX_mm, xc_mm);
                __MM yd_mm  = _MM_SUB(posY_mm, yc_mm);
                __MM xd2_mm = _MM_SQR(xd_mm);
                __MM yd2_mm = _MM_SQR(yd_mm);
                __MM xyd_mm = _MM_MUL(xd_mm,yd_mm);
        
                __MM det_mm = _MM_ADD(_MM_ADD(_MM_ADD(cp_mm,_MM_MUL(cx2_mm,xd2_mm)),
                                                            _MM_MUL(cy2_mm,yd2_mm)),
                                                            _MM_MUL(cxy_mm,xyd_mm));
                __MM flagPosDet = _MM_CMP(det_mm,zeros_mm,_CMP_GE_OS);
                if (_MM_TESTZ1(flagPosDet)) { continue; }

                __MM sqrtdet_mm = _MM_SQRT(_MM_AND(det_mm,flagPosDet));
                __MM cx1xcy1y_mm = _MM_ADD(_MM_ADD(_MM_MUL(cx1_mm,xd_mm),_MM_MUL(cy1_mm,yd_mm)),zc_mm);
                __MM tmpIdxMinZ_mm =  _MM_CEIL(_MM_SUB(cx1xcy1y_mm,sqrtdet_mm)); 
                __MM tmpIdxMaxZ_mm = _MM_FLOOR(_MM_ADD(cx1xcy1y_mm,sqrtdet_mm)); 
            
                __MM zstart_mm = _MM_MAX(tmpIdxMinZ_mm, zeros_mm); 
                __MM   zend_mm = _MM_MIN(tmpIdxMaxZ_mm,numzm1_mm); 

                __MM flagInRange = _MM_AND(flagPosDet,_MM_AND(
                        _MM_CMP(zstart_mm,numzm1_mm,_CMP_LE_OS), 
                        _MM_CMP(  zend_mm, zeros_mm,_CMP_GE_OS)));
                zstart_mm = _MM_BLENDV( numz_mm,zstart_mm,flagInRange);
                //zend_mm   = _MM_BLENDV(mones_mm,  zend_mm,flagInRange);
                zend_mm = _MM_BLENDV(zeros_mm,_MM_ADD(zend_mm,ones_mm),flagInRange);

                __MM zd_mm = _MM_SUB(zstart_mm,zc_mm);

                __MM xs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                                             _MM_MUL(yd_mm,Sinv_12_mm)),
                                             _MM_MUL(zd_mm,Sinv_13_mm));
                __MM ys_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                                             _MM_MUL(yd_mm,Sinv_22_mm)),
                                             _MM_MUL(zd_mm,Sinv_23_mm));
                __MM zs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                                             _MM_MUL(yd_mm,Sinv_23_mm)),
                                             _MM_MUL(zd_mm,Sinv_33_mm));
            
                __MM intensity0_mm = _MM_EXP(_MM_MUL(mhalf_mm,_MM_ADD(_MM_ADD(_MM_MUL(xd_mm,xs_mm),
                                                                              _MM_MUL(yd_mm,ys_mm)),
                                                                              _MM_MUL(zd_mm,zs_mm))));
                __MM az_mm = zs_mm;
                __MM fz0_mm = _MM_EXP(_MM_SUB(bz_mm,az_mm));

                //int tmpzmin,tmpzmax;
                //hmin(&tmpzmin,zstart_mm);
                //hmax(&tmpzmax,  zend_mm);
                int tmpzmin = hmin(zstart_mm);
                int tmpzmax = hmax(  zend_mm);
                for (int cz=tmpzmin; cz<tmpzmax; cz++) {
                    table_elipsoid_z[cz].push_back(cvalid);
                }             

                min_zstart_mm = _MM_MIN(min_zstart_mm,zstart_mm);
                max_zend_mm = _MM_MAX( max_zend_mm, zend_mm);
                idx_validelip[cvalid] = table_elipsoid_xy[tableidx][c];
                mem_int0_mm[cvalid] = intensity0_mm;
                mem_fz0_mm[cvalid] = fz0_mm;
                mem_int_mm[cvalid] = zeros_mm;
                mem_fz_mm[cvalid] = zeros_mm;
                idx_Zstart_mm[cvalid] = zstart_mm;
                idx_Zend_mm[cvalid] = zend_mm;
                cvalid++;           

            }
            unsigned int numvalid = cvalid;

            //int minZStart, maxZEnd;
            //hmin(&minZStart,min_zstart_mm);
            //hmax(  &maxZEnd,  max_zend_mm);
            int minZStart = hmin(min_zstart_mm);
            int   maxZEnd = hmax(  max_zend_mm);


            #ifdef MY_ASSERT
                t_init += __rdtsc() - t_start;
                t_start = __rdtsc();
            #endif


             /* loop for z */
            unsigned int cvxyz;
            for (int cz = minZStart; cz<maxZEnd; cz++ ){
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;

                if (table_elipsoid_z[cz].size()==0) { continue; }

                // update intensity of each elipsoid at the voxel
                __MM posZ_mm = _MM_SET1(cz);
                //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM flagstart_mm = _MM_CMP(idx_Zstart_mm[cvalid],posZ_mm,_CMP_EQ_OS);
                    __MM flagend_mm   = _MM_CMP(  idx_Zend_mm[cvalid],posZ_mm,_CMP_NEQ_OS);
                    __MM exp2bz_mm = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC + 17];
                
                    mem_int_mm[cvalid] = _MM_AND(_MM_ADD(_MM_MUL(mem_int_mm[cvalid],mem_fz_mm[cvalid]),
                                                         _MM_AND(mem_int0_mm[cvalid],flagstart_mm)),flagend_mm);
                    mem_fz_mm[cvalid]  =         _MM_ADD(_MM_MUL(mem_fz_mm[cvalid],exp2bz_mm),
                                                         _MM_AND(mem_fz0_mm[cvalid],flagstart_mm));
                
                }

                // when the voxel was masked, skip the following calculations                
                if (_MM_TESTZ1(mask_mm[cvxyz])) { continue; }
                
                // calc integrals 
                __MM intensity_mm = zeros_mm;
                //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    __MM pi_k_mm   = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  0];
                    __MM xc_mm     = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  1];
                    __MM yc_mm     = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  2];
                    __MM zc_mm     = params_proc_mm[idx_validelip[cvalid]*NUM_PARAMS_PROC +  3];
                    __MM tmpint_mm = _MM_MUL(mem_int_mm[cvalid],pi_k_mm);
                    intensity_mm = _MM_ADD(intensity_mm,tmpint_mm);
                
                    __MM tmpint2_mm = _MM_MUL(mem_int_mm[cvalid],tmpint_mm);
                    __MM xd_mm = _MM_SUB(posX_mm,xc_mm);
                    __MM yd_mm = _MM_SUB(posY_mm,yc_mm);
                    __MM zd_mm = _MM_SUB(posZ_mm,zc_mm);

                    new_S_11_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(xd_mm,xd_mm),tmpint2_mm);
                    new_S_12_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(xd_mm,yd_mm),tmpint2_mm);
                    new_S_13_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(xd_mm,zd_mm),tmpint2_mm);
                    new_S_22_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(yd_mm,yd_mm),tmpint2_mm);
                    new_S_23_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(yd_mm,zd_mm),tmpint2_mm);
                    new_S_33_tmp_mm[cvalid] = _MM_MUL(_MM_MUL(zd_mm,zd_mm),tmpint2_mm);                
                    sumrzg_tmp_mm  [cvalid] =                              tmpint2_mm ;

                }

                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                __MM residual_mm = _MM_SUB(im_orig_mm[cvxyz], intensity_mm);
                __MM rss_old_mm  = _MM_SQR(im_orig_mm[cvxyz]);
                __MM rss_tmp_mm  = _MM_SQR(residual_mm);
                twosub(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_old_mm);
                twosum(&rss_current_mm[idx_thread],&rss_current_err_mm[idx_thread],&rss_tmp_mm);

                __MM tmpratio_mm = _MM_DIV(residual_mm,intensity_mm);
                tmpratio_mm = _MM_AND(tmpratio_mm,_MM_CMP(intensity_mm,zeros_mm,_CMP_NEQ_OQ));
            
                 //for (int cvalid=0; cvalid<numvalid; cvalid++) {
                for (int ct=0; ct<table_elipsoid_z[cz].size(); ct++) {
                    cvalid = table_elipsoid_z[cz][ct];
                    unsigned int tmpidx = idx_validelip[cvalid] + omp_offset;
                    __MM tmp_S_11_tmp_mm   = _MM_MUL(new_S_11_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_S_12_tmp_mm   = _MM_MUL(new_S_12_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_S_13_tmp_mm   = _MM_MUL(new_S_13_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_S_22_tmp_mm   = _MM_MUL(new_S_22_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_S_23_tmp_mm   = _MM_MUL(new_S_23_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_S_33_tmp_mm   = _MM_MUL(new_S_33_tmp_mm[cvalid],tmpratio_mm);
                    __MM tmp_sumrzg_tmp_mm = _MM_MUL(  sumrzg_tmp_mm[cvalid],tmpratio_mm);
                    twosum(&new_S_11_mm[tmpidx],&new_S_11_err_mm[tmpidx],  &tmp_S_11_tmp_mm);
                    twosum(&new_S_12_mm[tmpidx],&new_S_12_err_mm[tmpidx],  &tmp_S_12_tmp_mm);
                    twosum(&new_S_13_mm[tmpidx],&new_S_13_err_mm[tmpidx],  &tmp_S_13_tmp_mm);
                    twosum(&new_S_22_mm[tmpidx],&new_S_22_err_mm[tmpidx],  &tmp_S_22_tmp_mm);
                    twosum(&new_S_23_mm[tmpidx],&new_S_23_err_mm[tmpidx],  &tmp_S_23_tmp_mm);
                    twosum(&new_S_33_mm[tmpidx],&new_S_33_err_mm[tmpidx],  &tmp_S_33_tmp_mm);
                    twosum(  &sumrzg_mm[tmpidx],  &sumrzg_err_mm[tmpidx],&tmp_sumrzg_tmp_mm);
                }
            }

            #ifdef MY_ASSERT
                t_loop += __rdtsc() - t_start;
            #endif

        }
    } // end of omp_for

    #ifdef _OPENMP

    // aggregate the local copies into global array
    #pragma omp for
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {
        for (unsigned int p=1; p<numthreads; p++) {

            twosum(&new_S_11_mm[celipsoid],&new_S_11_err_mm[celipsoid],&new_S_11_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_12_mm[celipsoid],&new_S_12_err_mm[celipsoid],&new_S_12_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_13_mm[celipsoid],&new_S_13_err_mm[celipsoid],&new_S_13_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_22_mm[celipsoid],&new_S_22_err_mm[celipsoid],&new_S_22_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_23_mm[celipsoid],&new_S_23_err_mm[celipsoid],&new_S_23_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_33_mm[celipsoid],&new_S_33_err_mm[celipsoid],&new_S_33_mm[numelipsoids*p+celipsoid]);
            twosum(  &sumrzg_mm[celipsoid],  &sumrzg_err_mm[celipsoid],  &sumrzg_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_11_mm[celipsoid],&new_S_11_err_mm[celipsoid],&new_S_11_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_12_mm[celipsoid],&new_S_12_err_mm[celipsoid],&new_S_12_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_13_mm[celipsoid],&new_S_13_err_mm[celipsoid],&new_S_13_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_22_mm[celipsoid],&new_S_22_err_mm[celipsoid],&new_S_22_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_23_mm[celipsoid],&new_S_23_err_mm[celipsoid],&new_S_23_err_mm[numelipsoids*p+celipsoid]);
            twosum(&new_S_33_mm[celipsoid],&new_S_33_err_mm[celipsoid],&new_S_33_err_mm[numelipsoids*p+celipsoid]);
            twosum(  &sumrzg_mm[celipsoid],  &sumrzg_err_mm[celipsoid],  &sumrzg_err_mm[numelipsoids*p+celipsoid]);
        }
    } // end of omp_for


    } // end of parallel

    // aggregate the local copies of rss_current_mm into global array
    for (unsigned int p=1; p<numthreads; p++) {
        twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_mm[p]);
        twosum(&rss_current_mm[0],&rss_current_err_mm[0],&rss_current_err_mm[p]);
    }
    #endif


    #ifdef MY_ASSERT
        t_start = __rdtsc();
    #endif    

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        /* sum up new_mu_* and set to tmpnum */
        SIMD_ELEMENT_TYPE tmpnum_S_11     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_12     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_13     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_22     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_23     = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_33     = 0.0;
        SIMD_ELEMENT_TYPE tmpden          = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_11_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_12_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_13_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_22_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_23_err = 0.0;
        SIMD_ELEMENT_TYPE tmpnum_S_33_err = 0.0;
        SIMD_ELEMENT_TYPE tmpden_err      = 0.0;
        for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
            twosum(&tmpnum_S_11,&tmpnum_S_11_err,&new_S_11[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_12,&tmpnum_S_12_err,&new_S_12[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_13,&tmpnum_S_13_err,&new_S_13[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_22,&tmpnum_S_22_err,&new_S_22[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_23,&tmpnum_S_23_err,&new_S_23[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_33,&tmpnum_S_33_err,&new_S_33[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpden,     &tmpden_err,       &sumrzg[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_11,&tmpnum_S_11_err,&new_S_11_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_12,&tmpnum_S_12_err,&new_S_12_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_13,&tmpnum_S_13_err,&new_S_13_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_22,&tmpnum_S_22_err,&new_S_22_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_23,&tmpnum_S_23_err,&new_S_23_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpnum_S_33,&tmpnum_S_33_err,&new_S_33_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
            twosum(&tmpden,     &tmpden_err,       &sumrzg_err[SIMD_VECTOR_LENGTH*celipsoid+counter]);
        }
        tmpnum_S_11 += tmpnum_S_11_err;
        tmpnum_S_12 += tmpnum_S_12_err;
        tmpnum_S_13 += tmpnum_S_13_err;
        tmpnum_S_22 += tmpnum_S_22_err;
        tmpnum_S_23 += tmpnum_S_23_err;
        tmpnum_S_33 += tmpnum_S_33_err;
        tmpden      += tmpden_err;

        /* update sigma */
        SIMD_ELEMENT_TYPE pi_k     = params[NUM_PARAMS_INPUT*celipsoid + 0];
        SIMD_ELEMENT_TYPE old_S_11 = params[NUM_PARAMS_INPUT*celipsoid + 4];
        SIMD_ELEMENT_TYPE old_S_12 = params[NUM_PARAMS_INPUT*celipsoid + 5];
        SIMD_ELEMENT_TYPE old_S_13 = params[NUM_PARAMS_INPUT*celipsoid + 6];
        SIMD_ELEMENT_TYPE old_S_22 = params[NUM_PARAMS_INPUT*celipsoid + 7];
        SIMD_ELEMENT_TYPE old_S_23 = params[NUM_PARAMS_INPUT*celipsoid + 8];
        SIMD_ELEMENT_TYPE old_S_33 = params[NUM_PARAMS_INPUT*celipsoid + 9];

        if (pi_k==0) { continue; }

        tmpnum_S_11 += old_S_11*pi_k*volume[celipsoid];
        tmpnum_S_12 += old_S_12*pi_k*volume[celipsoid];
        tmpnum_S_13 += old_S_13*pi_k*volume[celipsoid];
        tmpnum_S_22 += old_S_22*pi_k*volume[celipsoid];
        tmpnum_S_23 += old_S_23*pi_k*volume[celipsoid];
        tmpnum_S_33 += old_S_33*pi_k*volume[celipsoid];
        tmpden      +=          pi_k*volume[celipsoid];

        SIMD_ELEMENT_TYPE tmpden_inv = 1/tmpden;
        
        if (tmpnum_S_11*tmpden_inv<0 ||
            tmpnum_S_22*tmpden_inv<0 ||
            tmpnum_S_33*tmpden_inv<0 ) {
                continue;
        }

        SIMD_ELEMENT_TYPE detS = (    tmpnum_S_11*tmpnum_S_22*tmpnum_S_33 
                                  + 2*tmpnum_S_12*tmpnum_S_13*tmpnum_S_23
                                  -   tmpnum_S_11*tmpnum_S_23*tmpnum_S_23 
                                  -   tmpnum_S_12*tmpnum_S_12*tmpnum_S_33
                                  -   tmpnum_S_13*tmpnum_S_22*tmpnum_S_13)
                                 * tmpden_inv*tmpden_inv*tmpden_inv;
        if (detS<0) { continue; }
            


        params[NUM_PARAMS_INPUT*celipsoid + 4] = tmpnum_S_11 * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 5] = tmpnum_S_12 * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 6] = tmpnum_S_13 * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 7] = tmpnum_S_22 * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 8] = tmpnum_S_23 * tmpden_inv;
        params[NUM_PARAMS_INPUT*celipsoid + 9] = tmpnum_S_33 * tmpden_inv;
        
    }


    // update score
    SIMD_ELEMENT_TYPE tmpsum     = 0.0;
    SIMD_ELEMENT_TYPE tmpsum_err = 0.0;
    SIMD_ELEMENT_TYPE *rss_current     = (SIMD_ELEMENT_TYPE *)&rss_current_mm[0];
    SIMD_ELEMENT_TYPE *rss_current_err = (SIMD_ELEMENT_TYPE *)&rss_current_err_mm[0];
    for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
        twosum(&tmpsum,&tmpsum_err,&rss_current[counter]);
        twosum(&tmpsum,&tmpsum_err,&rss_current_err[counter]);
    }
    tmpsum += tmpsum_err;
    score = 1 + ((double)tmpsum) * rss_orig_inv;

    #ifdef MY_ASSERT
        t_update = __rdtsc() - t_start;
    #endif

    #ifdef MY_ASSERT    
        mexPrintf("em_optim.m2step.t_init: %I64d[clocks]\n", t_init);
        mexPrintf("em_optim.m2step.t_loop: %I64d[clocks]\n", t_loop);
        mexPrintf("em_optim.m2step.t_update: %I64d[clocks]\n", t_update);
    #endif

}


void calc_G() {
/* calc G (integral of gi*gj) for update pi_k */

    double *pG      = mxGetPr(prhs_mldivide[0]);
    


    // loop for elipsoid
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {        
        SIMD_ELEMENT_TYPE A1_11   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 0];
        SIMD_ELEMENT_TYPE A1_12   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 1];
        SIMD_ELEMENT_TYPE A1_13   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 2];
        SIMD_ELEMENT_TYPE A1_22   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 3];
        SIMD_ELEMENT_TYPE A1_23   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 4];
        SIMD_ELEMENT_TYPE A1_33   = params_calcG[celipsoid*NUM_PARAMS_CALCG + 5];
        SIMD_ELEMENT_TYPE B1_1    = params_calcG[celipsoid*NUM_PARAMS_CALCG + 6];
        SIMD_ELEMENT_TYPE B1_2    = params_calcG[celipsoid*NUM_PARAMS_CALCG + 7];
        SIMD_ELEMENT_TYPE B1_3    = params_calcG[celipsoid*NUM_PARAMS_CALCG + 8];
        SIMD_ELEMENT_TYPE C1      = params_calcG[celipsoid*NUM_PARAMS_CALCG + 9];


        for (unsigned int celipsoid2=0; celipsoid2<=celipsoid; celipsoid2++) {

            SIMD_ELEMENT_TYPE A2_11   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 0];
            SIMD_ELEMENT_TYPE A2_12   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 1];
            SIMD_ELEMENT_TYPE A2_13   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 2];
            SIMD_ELEMENT_TYPE A2_22   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 3];
            SIMD_ELEMENT_TYPE A2_23   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 4];
            SIMD_ELEMENT_TYPE A2_33   = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 5];
            SIMD_ELEMENT_TYPE B2_1    = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 6];
            SIMD_ELEMENT_TYPE B2_2    = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 7];
            SIMD_ELEMENT_TYPE B2_3    = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 8];
            SIMD_ELEMENT_TYPE C2      = params_calcG[celipsoid2*NUM_PARAMS_CALCG + 9];

            if (   idx_VXstart[celipsoid ] > idx_VXend[celipsoid2]
                || idx_VXstart[celipsoid2] > idx_VXend[celipsoid ] 
                ||  idx_Ystart[celipsoid ] >  idx_Yend[celipsoid2]
                ||  idx_Ystart[celipsoid2] >  idx_Yend[celipsoid ] ) { continue; }


            SIMD_ELEMENT_TYPE A_11 = A1_11 + A2_11;
            SIMD_ELEMENT_TYPE A_22 = A1_22 + A2_22;
            SIMD_ELEMENT_TYPE A_33 = A1_33 + A2_33;
            SIMD_ELEMENT_TYPE A_12 = A1_12 + A2_12;
            SIMD_ELEMENT_TYPE A_13 = A1_13 + A2_13;
            SIMD_ELEMENT_TYPE A_23 = A1_23 + A2_23;

            SIMD_ELEMENT_TYPE B_1 = B1_1 + B2_1;
            SIMD_ELEMENT_TYPE B_2 = B1_2 + B2_2;
            SIMD_ELEMENT_TYPE B_3 = B1_3 + B2_3;

            SIMD_ELEMENT_TYPE C = C1 + C2;

            
            SIMD_ELEMENT_TYPE detA =    A_11*A_22*A_33 
                                    + 2*A_12*A_13*A_23
                                    -   A_11*A_23*A_23 
                                    -   A_12*A_12*A_33
                                    -   A_13*A_22*A_13;

            SIMD_ELEMENT_TYPE detA_inv = 1/detA;

            SIMD_ELEMENT_TYPE Ainv_11 = A_22*A_33 - A_23*A_23;
            SIMD_ELEMENT_TYPE Ainv_22 = A_11*A_33 - A_13*A_13;
            SIMD_ELEMENT_TYPE Ainv_33 = A_11*A_22 - A_12*A_12;
            SIMD_ELEMENT_TYPE Ainv_12 = A_13*A_23 - A_12*A_33;
            SIMD_ELEMENT_TYPE Ainv_13 = A_12*A_23 - A_13*A_22;
            SIMD_ELEMENT_TYPE Ainv_23 = A_12*A_13 - A_11*A_23;

            
            SIMD_ELEMENT_TYPE BAB = (  ( B_1*Ainv_11 + B_2*Ainv_12 + B_3*Ainv_13 ) * B_1
                                     + ( B_1*Ainv_12 + B_2*Ainv_22 + B_3*Ainv_23 ) * B_2
                                     + ( B_1*Ainv_13 + B_2*Ainv_23 + B_3*Ainv_33 ) * B_3 ) 
                                    * detA_inv;


            SIMD_ELEMENT_TYPE gg = sqrt(8*M_PI*M_PI*M_PI*detA_inv)*exp(0.5*(BAB-C));

            pG[celipsoid  * numelipsoids + celipsoid2] = (double) gg;
            pG[celipsoid2 * numelipsoids + celipsoid ] = (double) gg;

        }
    }
}


void calc_H() {
    /* calc H (integral of gi*im) for update pi_k */

    #ifdef MY_ASSERT    
        unsigned __int64 t_start, t_init=0, t_loop=0, t_update=0;
    #endif

    double *pH = mxGetPr(prhs_mldivide[1]);
    
    // loop for elipsoid
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        unsigned int offset = celipsoid*NUM_PARAMS_PROC;

        __MM xc_mm      = params_proc_mm[offset +  1];
        __MM yc_mm      = params_proc_mm[offset +  2];
        __MM zc_mm      = params_proc_mm[offset +  3];
        __MM Sinv_11_mm = params_proc_mm[offset +  4];
        __MM Sinv_12_mm = params_proc_mm[offset +  5];
        __MM Sinv_13_mm = params_proc_mm[offset +  6];
        __MM Sinv_22_mm = params_proc_mm[offset +  7];
        __MM Sinv_23_mm = params_proc_mm[offset +  8];
        __MM Sinv_33_mm = params_proc_mm[offset +  9];
        __MM cx1_mm     = params_proc_mm[offset + 10];
        __MM cy1_mm     = params_proc_mm[offset + 11];
        __MM cx2_mm     = params_proc_mm[offset + 12];
        __MM cy2_mm     = params_proc_mm[offset + 13];
        __MM cxy_mm     = params_proc_mm[offset + 14];
        __MM cp_mm      = params_proc_mm[offset + 15];
        __MM bz_mm      = params_proc_mm[offset + 16];
        __MM exp2bz_mm  = params_proc_mm[offset + 17];

        __MM hh_mm     = zeros_mm;
        __MM hh_err_mm = zeros_mm;

        

        for (unsigned int cy = idx_Ystart[celipsoid]; cy<=idx_Yend[celipsoid]; cy++) {
            __MM yd_mm  = _MM_SUB(_MM_SET1(cy), yc_mm);
            __MM yd2_mm = _MM_SQR(yd_mm);


            for (unsigned int cvx = idx_VXstart[celipsoid]; cvx<=idx_VXend[celipsoid]; cvx++) {
                
                #ifdef MY_ASSERT    
                    t_start = __rdtsc();
                #endif


                __MM xd_mm = _MM_SUB(_MM_ADD(_MM_SET1(cvx*SIMD_VECTOR_LENGTH),addx_mm), xc_mm);
                __MM xd2_mm = _MM_SQR(xd_mm);
                __MM xyd_mm = _MM_MUL(xd_mm,yd_mm);

                __MM det_mm = _MM_ADD(_MM_ADD(_MM_ADD(cp_mm,_MM_MUL(cx2_mm,xd2_mm)),
                                                            _MM_MUL(cy2_mm,yd2_mm)),
                                                            _MM_MUL(cxy_mm,xyd_mm));
                __MM flagPosDet = _MM_CMP(det_mm,zeros_mm,_CMP_GE_OS);
                if (_MM_TESTZ1(flagPosDet)) { continue; }

                __MM sqrtdet_mm = _MM_SQRT(_MM_AND(det_mm,flagPosDet));
                __MM cx1xcy1y_mm = _MM_ADD(_MM_ADD(_MM_MUL(cx1_mm,xd_mm),_MM_MUL(cy1_mm,yd_mm)),zc_mm);
                __MM tmpIdxMinZ_mm =  _MM_CEIL(_MM_SUB(cx1xcy1y_mm,sqrtdet_mm)); 
                __MM tmpIdxMaxZ_mm = _MM_FLOOR(_MM_ADD(cx1xcy1y_mm,sqrtdet_mm)); 
            
                __MM zstart_mm = _MM_MAX(tmpIdxMinZ_mm, zeros_mm); 
                __MM   zend_mm = _MM_MIN(tmpIdxMaxZ_mm,numzm1_mm); 

                __MM flagInRange = _MM_AND(flagPosDet,_MM_AND(
                        _MM_CMP(zstart_mm,numzm1_mm,_CMP_LE_OS), 
                        _MM_CMP(  zend_mm, zeros_mm,_CMP_GE_OS)));
                zstart_mm = _MM_BLENDV( numz_mm,zstart_mm,flagInRange);
                //zend_mm   = _MM_BLENDV(mones_mm,  zend_mm,flagInRange);
                zend_mm = _MM_BLENDV(zeros_mm,_MM_ADD(zend_mm,ones_mm),flagInRange);

                __MM zd_mm = _MM_SUB(zstart_mm,zc_mm);

                __MM xs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_11_mm),
                                             _MM_MUL(yd_mm,Sinv_12_mm)),
                                             _MM_MUL(zd_mm,Sinv_13_mm));
                __MM ys_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_12_mm),
                                             _MM_MUL(yd_mm,Sinv_22_mm)),
                                             _MM_MUL(zd_mm,Sinv_23_mm));
                __MM zs_mm = _MM_ADD(_MM_ADD(_MM_MUL(xd_mm,Sinv_13_mm),
                                             _MM_MUL(yd_mm,Sinv_23_mm)),
                                             _MM_MUL(zd_mm,Sinv_33_mm));
            
                __MM intensity0_mm = _MM_EXP(_MM_MUL(mhalf_mm,_MM_ADD(_MM_ADD(_MM_MUL(xd_mm,xs_mm),
                                                                              _MM_MUL(yd_mm,ys_mm)),
                                                                              _MM_MUL(zd_mm,zs_mm))));
                __MM az_mm = zs_mm;
                __MM fz0_mm = _MM_EXP(_MM_SUB(bz_mm,az_mm));


                int minZStart = hmin(zstart_mm);
                int   maxZEnd = hmax(  zend_mm);


                //__MMI zstart_mmi =_MM_CVTF2I(zstart_mm);
                //__MMI   zend_mmi =_MM_CVTF2I(  zend_mm);
                //int *zstart = (int *)&zstart_mmi;
                //int   *zend = (int *)&zend_mmi;

                //#if defined(ENABLE_SINGLE_PRECISION) // horizontal min and max for single float mm256
                //    int minZStart = __min(__min(__min(zstart[0],zstart[1]),
                //                                __min(zstart[2],zstart[3])),
                //                          __min(__min(zstart[4],zstart[5]),
                //                                __min(zstart[6],zstart[7])));
                //    int maxZEnd = __max(__max(__max(zend[0],zend[1]),
                //                              __max(zend[2],zend[3])),
                //                        __max(__max(zend[4],zend[5]),
                //                              __max(zend[6],zend[7])));
                //#else // horizontal min and max for double float mm256d
                //    int minZStart = __min(__min(zstart[0],zstart[1]),
                //                          __min(zstart[2],zstart[3]));                                            
                //    int maxZEnd = __max(__max(zend[0],zend[1]),
                //                        __max(zend[2],zend[3]));
                //#endif

                #ifdef MY_ASSERT
                    t_init += __rdtsc() - t_start;
                    t_start = __rdtsc();
                #endif

                // /* set the initial values to the tables */
                //SIMD_ELEMENT_TYPE *intensity0 = (SIMD_ELEMENT_TYPE *)&intensity0_mm;
                //SIMD_ELEMENT_TYPE *fz0        = (SIMD_ELEMENT_TYPE *)&fz0_mm;
                //for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
                //    im_init_int[ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = intensity0[counter];
                //    im_init_fz [ zstart[counter] *SIMD_VECTOR_LENGTH+counter] =        fz0[counter];
                //    im_end_int [(1+zend[counter])*SIMD_VECTOR_LENGTH+counter] = mask_false;                    
                //}


                 /* loop for z */
                __MM intensity_mm = zeros_mm;
                __MM fz_mm = zeros_mm;
                unsigned int cvxyz;
                //for (int cz = minZStart; cz<=maxZEnd; cz++ ){
                for (int cz = minZStart; cz<maxZEnd; cz++ ){
                    /* z方向の更新式 */
                    //intensity_mm = _MM_AND(_MM_ADD(_MM_MUL(intensity_mm,fz_mm),imInitInt_mm[cz]),imEndInt_mm[cz]);
                    //fz_mm        =         _MM_ADD(_MM_MUL(fz_mm,exp2bz_mm),   imInitFz_mm[cz]);
                    __MM posZ_mm = _MM_SET1(cz); 
                    __MM flagstart_mm = _MM_CMP(zstart_mm,posZ_mm,_CMP_EQ_OS);
                    __MM flagend_mm   = _MM_CMP(  zend_mm,posZ_mm,_CMP_NEQ_OS);
                
                    intensity_mm = _MM_AND(_MM_ADD(_MM_MUL(intensity_mm,fz_mm),_MM_AND(intensity0_mm,flagstart_mm)),flagend_mm);
                    fz_mm        =         _MM_ADD(_MM_MUL(fz_mm,exp2bz_mm),   _MM_AND(fz0_mm,flagstart_mm));
                
                    cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                    __MM hh_tmp_mm = _MM_AND(_MM_MUL(im_orig_mm[cvxyz],intensity_mm),mask_mm[cvxyz]);
                    twosum(&hh_mm,&hh_err_mm,&hh_tmp_mm);
                }

                ///* reset initial value tables */
                //for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
                //    im_init_int[ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = 0.0;
                //    im_init_fz [ zstart[counter] *SIMD_VECTOR_LENGTH+counter] = 0.0;
                //    im_end_int [(1+zend[counter])*SIMD_VECTOR_LENGTH+counter] = mask_true;
                //}

                #ifdef MY_ASSERT    
                    t_loop += __rdtsc() - t_start;
                #endif

            }
        }

        #ifdef MY_ASSERT
            t_start = __rdtsc();
        #endif

        SIMD_ELEMENT_TYPE *hh     = (SIMD_ELEMENT_TYPE *)&hh_mm;
        SIMD_ELEMENT_TYPE *hh_err = (SIMD_ELEMENT_TYPE *)&hh_err_mm;
        SIMD_ELEMENT_TYPE tmpsum = 0.0;
        SIMD_ELEMENT_TYPE tmperr = 0.0;
        for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
            twosum(&tmpsum,&tmperr,&hh[counter]);
            twosum(&tmpsum,&tmperr,&hh_err[counter]);
        }
        pH[celipsoid] = (double) (tmpsum+tmperr);

        #ifdef MY_ASSERT
            t_update = __rdtsc() - t_start;
        #endif

    }

    #ifdef MY_ASSERT    
        mexPrintf("em_optim.estep.calc_H.t_init: %I64d[clocks]\n", t_init);
        mexPrintf("em_optim.estep.calc_H.t_loop: %I64d[clocks]\n", t_loop);
        mexPrintf("em_optim.estep.calc_H.t_update: %I64d[clocks]\n", t_update);
    #endif

}


void setparam() {    

    SIMD_ELEMENT_TYPE c = __min(2*thrdist*thrdist,-FLT_EXP_MIN_INPUT);

    for (unsigned int cvxy=0; cvxy<table_elipsoid_xy.size(); cvxy++) {
            table_elipsoid_xy[cvxy].clear();
    }

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        SIMD_ELEMENT_TYPE pi   = params[NUM_PARAMS_INPUT*celipsoid + 0];
        SIMD_ELEMENT_TYPE xc   = params[NUM_PARAMS_INPUT*celipsoid + 1];
        SIMD_ELEMENT_TYPE yc   = params[NUM_PARAMS_INPUT*celipsoid + 2];
        SIMD_ELEMENT_TYPE zc   = params[NUM_PARAMS_INPUT*celipsoid + 3];
        SIMD_ELEMENT_TYPE S_11 = params[NUM_PARAMS_INPUT*celipsoid + 4];
        SIMD_ELEMENT_TYPE S_12 = params[NUM_PARAMS_INPUT*celipsoid + 5];
        SIMD_ELEMENT_TYPE S_13 = params[NUM_PARAMS_INPUT*celipsoid + 6];
        SIMD_ELEMENT_TYPE S_22 = params[NUM_PARAMS_INPUT*celipsoid + 7];
        SIMD_ELEMENT_TYPE S_23 = params[NUM_PARAMS_INPUT*celipsoid + 8];
        SIMD_ELEMENT_TYPE S_33 = params[NUM_PARAMS_INPUT*celipsoid + 9];


        SIMD_ELEMENT_TYPE detS =    S_11*S_22*S_33 
                                + 2*S_12*S_13*S_23
                                -   S_11*S_23*S_23 
                                -   S_12*S_12*S_33
                                -   S_13*S_22*S_13;

        SIMD_ELEMENT_TYPE detS_inv = 1/detS;

        SIMD_ELEMENT_TYPE Sinv_11 = detS_inv*(S_22*S_33 - S_23*S_23);
        SIMD_ELEMENT_TYPE Sinv_22 = detS_inv*(S_11*S_33 - S_13*S_13);
        SIMD_ELEMENT_TYPE Sinv_33 = detS_inv*(S_11*S_22 - S_12*S_12);
        SIMD_ELEMENT_TYPE Sinv_12 = detS_inv*(S_13*S_23 - S_12*S_33);
        SIMD_ELEMENT_TYPE Sinv_13 = detS_inv*(S_12*S_23 - S_13*S_22);
        SIMD_ELEMENT_TYPE Sinv_23 = detS_inv*(S_12*S_13 - S_11*S_23);


        SIMD_ELEMENT_TYPE B_1 = xc*Sinv_11 + yc*Sinv_12 + zc*Sinv_13;
        SIMD_ELEMENT_TYPE B_2 = xc*Sinv_12 + yc*Sinv_22 + zc*Sinv_23;
        SIMD_ELEMENT_TYPE B_3 = xc*Sinv_13 + yc*Sinv_23 + zc*Sinv_33;
        SIMD_ELEMENT_TYPE C_1 = B_1*xc + B_2*yc + B_3*zc;
        

        SIMD_ELEMENT_TYPE exp2bz = exp(-Sinv_33);

        SIMD_ELEMENT_TYPE Sinv_33_inv = 1/Sinv_33;
        SIMD_ELEMENT_TYPE Sp_11 = Sinv_11 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_12 = Sinv_12 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_13 = Sinv_13 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_22 = Sinv_22 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_23 = Sinv_23 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE Sp_33 = Sinv_33 * Sinv_33_inv;
        SIMD_ELEMENT_TYPE cp    = c       * Sinv_33_inv;

        SIMD_ELEMENT_TYPE cx2 = Sp_13*Sp_13 - Sp_11;
        SIMD_ELEMENT_TYPE cy2 = Sp_23*Sp_23 - Sp_22;
        SIMD_ELEMENT_TYPE cxy = Sp_13*Sp_23 - Sp_12;
        SIMD_ELEMENT_TYPE comden_inv = cp / ( cxy*cxy - cx2*cy2 );

        SIMD_ELEMENT_TYPE xsqrt = sqrt( cy2 * comden_inv );
        SIMD_ELEMENT_TYPE ysqrt = sqrt( cx2 * comden_inv );
        
        unsigned int vxstart = (unsigned int) (__min(__max( ceil(xc - xsqrt),0),numx-1)/SIMD_VECTOR_LENGTH);
        unsigned int vxend   = (unsigned int) (__min(__max(floor(xc + xsqrt),0),numx-1)/SIMD_VECTOR_LENGTH);
        unsigned int ystart  = (unsigned int) (__min(__max( ceil(yc - ysqrt),0),numy-1));        
        unsigned int yend    = (unsigned int) (__min(__max(floor(yc + ysqrt),0),numy-1));


        for (unsigned int cy=ystart; cy<=yend; cy++) {
            for (unsigned int cvx=vxstart; cvx<=vxend; cvx++) {
                table_elipsoid_xy[cy*numVectorX+cvx].push_back(celipsoid);
            }
        }


        volume[celipsoid] = sqrt(8*M_PI*M_PI*M_PI)*sqrt(detS);

        idx_VXstart[celipsoid] = vxstart;
        idx_VXend  [celipsoid] = vxend;
        idx_Ystart [celipsoid] = ystart;
        idx_Yend   [celipsoid] = yend;

        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  0] = _MM_SET1(pi);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  1] = _MM_SET1(xc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  2] = _MM_SET1(yc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  3] = _MM_SET1(zc);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  4] = _MM_SET1(Sinv_11);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  5] = _MM_SET1(Sinv_12);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  6] = _MM_SET1(Sinv_13);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  7] = _MM_SET1(Sinv_22); 
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  8] = _MM_SET1(Sinv_23);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid +  9] = _MM_SET1(Sinv_33);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 10] = _MM_SET1(-Sp_13); // cx1
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 11] = _MM_SET1(-Sp_23); // cy1
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 12] = _MM_SET1(cx2); 
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 13] = _MM_SET1(cy2);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 14] = _MM_SET1(2*cxy);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 15] = _MM_SET1(cp);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 16] = _MM_SET1(-0.5*Sinv_33); // -bz
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 17] = _MM_SET1(exp2bz);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 18] = _MM_SET1(S_11);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 19] = _MM_SET1(S_12);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 20] = _MM_SET1(S_13);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 21] = _MM_SET1(S_22);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 22] = _MM_SET1(S_23);
        params_proc_mm[NUM_PARAMS_PROC*celipsoid + 23] = _MM_SET1(S_33);

        params_calcG[NUM_PARAMS_CALCG*celipsoid + 0] = Sinv_11; // A_11
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 1] = Sinv_12; // A_12
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 2] = Sinv_13; // A_13
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 3] = Sinv_22; // A_22
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 4] = Sinv_23; // A_23
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 5] = Sinv_33; // A_33
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 6] = B_1; // B_1
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 7] = B_2; // B_2
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 8] = B_3; // B_3
        params_calcG[NUM_PARAMS_CALCG*celipsoid + 9] = C_1; // C
        

    }

}


void setparam_init(const mxArray *prhs0){

	/* get the number of parameters and parameter sets */
	unsigned int m = (unsigned int)mxGetM(prhs0);
	unsigned int n = (unsigned int)mxGetN(prhs0);

    int maxnumthreads = 1;
    #ifdef _OPENMP
        maxnumthreads = omp_get_max_threads();
    #endif

	if ( (m*n) != NUM_PARAMS_INPUT*numelipsoids) {

        numelipsoids = (m*n) / NUM_PARAMS_INPUT;
        
        free_params(); // free memory for params 

        idx_validelip = (int *)_aligned_malloc(sizeof(int)*numelipsoids*maxnumthreads,MEMORY_ALIGNMENT);
        idx_minZstart = (int *)_aligned_malloc(sizeof(int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_maxZend   = (int *)_aligned_malloc(sizeof(int)*numelipsoids,MEMORY_ALIGNMENT);


        idx_VXstart = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_Ystart  = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_VXend   = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        idx_Yend    = (unsigned int *)_aligned_malloc(sizeof(unsigned int)*numelipsoids,MEMORY_ALIGNMENT);
        
        size_t tmpsize = sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*numelipsoids*maxnumthreads;
        idx_Zstart   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        idx_Zend     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_int0     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_fz0      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_int      = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        mem_fz       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33_err = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg_err   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_x_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_y_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_mu_z_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_11_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_12_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_13_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_22_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_23_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        new_S_33_tmp = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        sumrzg_tmp   = (SIMD_ELEMENT_TYPE *)_aligned_malloc(tmpsize,MEMORY_ALIGNMENT);
        
        params_proc  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*numelipsoids*NUM_PARAMS_PROC,MEMORY_ALIGNMENT);
        params       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_INPUT,MEMORY_ALIGNMENT);
        volume       = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids,MEMORY_ALIGNMENT);
        params_calcG = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numelipsoids*NUM_PARAMS_CALCG,MEMORY_ALIGNMENT);
     

        /* SIMD配列へ変換 */
        idx_Zstart_mm = (__MM *) idx_Zstart;
        idx_Zend_mm = (__MM *) idx_Zend;
        mem_int0_mm = (__MM *) mem_int0;
        mem_fz0_mm  = (__MM *) mem_fz0;
        mem_int_mm = (__MM *) mem_int;
        mem_fz_mm  = (__MM *) mem_fz;
        new_mu_x_mm = (__MM *) new_mu_x;
        new_mu_y_mm = (__MM *) new_mu_y;
        new_mu_z_mm = (__MM *) new_mu_z;
        new_S_11_mm = (__MM *) new_S_11;
        new_S_12_mm = (__MM *) new_S_12;
        new_S_13_mm = (__MM *) new_S_13;
        new_S_22_mm = (__MM *) new_S_22;
        new_S_23_mm = (__MM *) new_S_23;
        new_S_33_mm = (__MM *) new_S_33;
        sumrzg_mm = (__MM *) sumrzg;
        new_mu_x_err_mm = (__MM *) new_mu_x_err;
        new_mu_y_err_mm = (__MM *) new_mu_y_err;
        new_mu_z_err_mm = (__MM *) new_mu_z_err;
        new_S_11_err_mm = (__MM *) new_S_11_err;
        new_S_12_err_mm = (__MM *) new_S_12_err;
        new_S_13_err_mm = (__MM *) new_S_13_err;
        new_S_22_err_mm = (__MM *) new_S_22_err;
        new_S_23_err_mm = (__MM *) new_S_23_err;
        new_S_33_err_mm = (__MM *) new_S_33_err;
        sumrzg_err_mm = (__MM *) sumrzg_err;
        new_mu_x_tmp_mm = (__MM *) new_mu_x_tmp;
        new_mu_y_tmp_mm = (__MM *) new_mu_y_tmp;
        new_mu_z_tmp_mm = (__MM *) new_mu_z_tmp;
        new_S_11_tmp_mm = (__MM *) new_S_11_tmp;
        new_S_12_tmp_mm = (__MM *) new_S_12_tmp;
        new_S_13_tmp_mm = (__MM *) new_S_13_tmp;
        new_S_22_tmp_mm = (__MM *) new_S_22_tmp;
        new_S_23_tmp_mm = (__MM *) new_S_23_tmp;
        new_S_33_tmp_mm = (__MM *) new_S_33_tmp;
        sumrzg_tmp_mm = (__MM *) sumrzg_tmp;
        params_proc_mm = (__MM *) params_proc;

	}

	double *pin = mxGetPr(prhs0);

    for (unsigned int celipsoid=0; celipsoid<numelipsoids; celipsoid++) {

        params[NUM_PARAMS_INPUT*celipsoid + 0] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 0]; // pi
        params[NUM_PARAMS_INPUT*celipsoid + 1] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 1] - 1; // xc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 2] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 2] - 1; // yc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 3] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 3] - 1; // zc, 0-start
        params[NUM_PARAMS_INPUT*celipsoid + 4] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 4]; // S_11
        params[NUM_PARAMS_INPUT*celipsoid + 5] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 5]; // S_12
        params[NUM_PARAMS_INPUT*celipsoid + 6] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 6]; // S_13
        params[NUM_PARAMS_INPUT*celipsoid + 7] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 7]; // S_22
        params[NUM_PARAMS_INPUT*celipsoid + 8] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 8]; // S_23
        params[NUM_PARAMS_INPUT*celipsoid + 9] = (SIMD_ELEMENT_TYPE) pin[NUM_PARAMS_INPUT*celipsoid + 9]; // S_33

    }

}


void free_im(){

    /* グローバル変数用のメモリを解放 */
    
    _aligned_free(mask);
    _aligned_free(im_orig);
    _aligned_free(im_synth);
    //_aligned_free(im_init_int);
    //_aligned_free(im_init_fz);
    //_aligned_free(im_end_int);
    //_aligned_free(flag_invalid_allZ);
    _aligned_free(rss_current_mm);
    _aligned_free(rss_current_err_mm);
    
	/* 念のためグローバル変数用のポインタにnullを代入 */

    mask = NULL;
    im_orig = NULL;
    im_synth = NULL;
    //im_init_int = NULL;
    //im_init_fz = NULL;
    //im_end_int = NULL;
    //flag_invalid_allZ = NULL;
    
    mask_mm     = NULL;
    im_orig_mm  = NULL;
    im_synth_mm  = NULL;
    //imInitInt_mm = NULL;
    //imInitFz_mm  = NULL;
    //imEndInt_mm  = NULL;
    rss_current_mm = NULL;
    rss_current_err_mm = NULL;
   
    /* detect memory leaks using msvc */
    #ifdef _CRTDBG_MAP_ALLOC
    _CrtDumpMemoryLeaks();
    #endif

}


void free_params() {

	/* グローバル変数用のメモリを解放 */
    _aligned_free(idx_validelip);
    _aligned_free(idx_minZstart);
    _aligned_free(idx_maxZend);

    _aligned_free(idx_VXstart);
    _aligned_free(idx_Ystart);
    _aligned_free(idx_VXend);
    _aligned_free(idx_Yend);
    
    _aligned_free(idx_Zstart);
    _aligned_free(idx_Zend);
    _aligned_free(mem_int0);
    _aligned_free(mem_fz0);
    _aligned_free(mem_int);
    _aligned_free(mem_fz);
    _aligned_free(new_mu_x);
    _aligned_free(new_mu_y);
    _aligned_free(new_mu_z);
    _aligned_free(new_S_11);
    _aligned_free(new_S_12);
    _aligned_free(new_S_13);
    _aligned_free(new_S_22);
    _aligned_free(new_S_23);
    _aligned_free(new_S_33);
    _aligned_free(sumrzg);
    _aligned_free(new_mu_x_err);
    _aligned_free(new_mu_y_err);
    _aligned_free(new_mu_z_err);
    _aligned_free(new_S_11_err);
    _aligned_free(new_S_12_err);
    _aligned_free(new_S_13_err);
    _aligned_free(new_S_22_err);
    _aligned_free(new_S_23_err);
    _aligned_free(new_S_33_err);
    _aligned_free(sumrzg_err);
    _aligned_free(new_mu_x_tmp);
    _aligned_free(new_mu_y_tmp);
    _aligned_free(new_mu_z_tmp);
    _aligned_free(new_S_11_tmp);
    _aligned_free(new_S_12_tmp);
    _aligned_free(new_S_13_tmp);
    _aligned_free(new_S_22_tmp);
    _aligned_free(new_S_23_tmp);
    _aligned_free(new_S_33_tmp);
    _aligned_free(sumrzg_tmp);
    _aligned_free(params_proc);
    _aligned_free(params);
    _aligned_free(volume);
    _aligned_free(params_calcG);

	/* 念のためグローバル変数用のポインタにnullを代入 */
    idx_validelip = NULL;
    idx_minZstart = NULL;
    idx_maxZend = NULL;

    idx_VXstart = NULL;
    idx_Ystart = NULL;
    idx_VXend = NULL;
    idx_Yend = NULL;
    
    idx_Zstart = NULL;
    idx_Zend = NULL;
    mem_int0 = NULL;
    mem_fz0 = NULL;
    mem_int = NULL;
    mem_fz = NULL;
    new_mu_x = NULL;
    new_mu_y = NULL;
    new_mu_z = NULL;
    new_S_11 = NULL;
    new_S_12 = NULL;
    new_S_13 = NULL;
    new_S_22 = NULL;
    new_S_23 = NULL;
    new_S_33 = NULL;
    sumrzg = NULL;
    new_mu_x_err = NULL;
    new_mu_y_err = NULL;
    new_mu_z_err = NULL;
    new_S_11_err = NULL;
    new_S_12_err = NULL;
    new_S_13_err = NULL;
    new_S_22_err = NULL;
    new_S_23_err = NULL;
    new_S_33_err = NULL;
    sumrzg_err = NULL;
    new_mu_x_tmp = NULL;
    new_mu_y_tmp = NULL;
    new_mu_z_tmp = NULL;
    new_S_11_tmp = NULL;
    new_S_12_tmp = NULL;
    new_S_13_tmp = NULL;
    new_S_22_tmp = NULL;
    new_S_23_tmp = NULL;
    new_S_33_tmp = NULL;
    sumrzg_tmp = NULL;
    params_proc = NULL;
    params = NULL;
    volume = NULL;
    params_calcG = NULL;

    /* SIMD配列のポインタにもnullを代入 */
    idx_Zstart_mm = NULL;
    idx_Zend_mm = NULL;
    mem_int0_mm = NULL;
    mem_fz0_mm = NULL;
    mem_int_mm = NULL;
    mem_fz_mm = NULL;
    new_mu_x_mm = NULL;
    new_mu_y_mm = NULL;
    new_mu_z_mm = NULL;
    new_S_11_mm = NULL;
    new_S_12_mm = NULL;
    new_S_13_mm = NULL;
    new_S_22_mm = NULL;
    new_S_23_mm = NULL;
    new_S_33_mm = NULL;
    sumrzg_mm = NULL;
    new_mu_x_err_mm = NULL;
    new_mu_y_err_mm = NULL;
    new_mu_z_err_mm = NULL;
    new_S_11_err_mm = NULL;
    new_S_12_err_mm = NULL;
    new_S_13_err_mm = NULL;
    new_S_22_err_mm = NULL;
    new_S_23_err_mm = NULL;
    new_S_33_err_mm = NULL;
    sumrzg_err_mm = NULL;
    new_mu_x_tmp_mm = NULL;
    new_mu_y_tmp_mm = NULL;
    new_mu_z_tmp_mm = NULL;
    new_S_11_tmp_mm = NULL;
    new_S_12_tmp_mm = NULL;
    new_S_13_tmp_mm = NULL;
    new_S_22_tmp_mm = NULL;
    new_S_23_tmp_mm = NULL;
    new_S_33_tmp_mm = NULL;
    sumrzg_tmp_mm = NULL;
    params_proc_mm = NULL;
    
    /* detect memory leaks using msvc */
    #ifdef _CRTDBG_MAP_ALLOC
    _CrtDumpMemoryLeaks();
    #endif

}


static void closefun(void) {
    free_im();
    free_params();
}


void initialize(int nrhs, const mxArray *prhs[]) {
	/* 変数宣言 */
	
    /* mex終了時に呼び出してメモリをfree()する関数を登録 */
	mexAtExit(closefun);
	
    /* 既に確保している場合にメモリを開放 */
	free_im();		

    /* processing option settings */

    // set default values
    thrint = -FLT_MAX_PRECISION;
    thrdist = sqrt(-FLT_EXP_MIN_INPUT);
    tol = 5e-5;
    maxiter = 1000;
    fixmu = 0;    

    // set the specified values
    if (nrhs==3 && mxIsStruct(prhs[2])) {
        mxArray *tmparr;
        if (mxGetFieldNumber(prhs[2],"thrint")>=0) {
            tmparr = mxGetField(prhs[2],0,"thrint");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                thrint = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"thrdist")>=0) {
            tmparr = mxGetField(prhs[2],0,"thrdist");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                thrdist = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"tol")>=0) {
            tmparr = mxGetField(prhs[2],0,"tol");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                tol = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"maxiter")>=0) {
            tmparr = mxGetField(prhs[2],0,"maxiter");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                maxiter = mxGetScalar(tmparr);
            }
        }
        if (mxGetFieldNumber(prhs[2],"fixmu")>=0) {
            tmparr = mxGetField(prhs[2],0,"fixmu");
            if (mxGetM(tmparr)*mxGetN(tmparr)==1) {
                fixmu = mxGetScalar(tmparr);
            }
        }
    }
	
	/* get the size of the image: [numx,numy,numz] = size(im)*/
	dims = mxGetDimensions(prhs[1]);
	numdims = mxGetNumberOfDimensions(prhs[1]);
	numx = (unsigned int)dims[0]; // size(im,1)
	numy = (unsigned int)dims[1]; // size(im,2)
	if ((unsigned int) numdims >= 3) {
		numz = (unsigned int)dims[2]; // size(im,3)
	} else {
		numz = 1;  // size(im,3)
	}
    numz_mm = _MM_SET1(numz);
    numzm1_mm = _MM_SET1(numz-1);
	
    padlength = (SIMD_VECTOR_LENGTH - (numx % SIMD_VECTOR_LENGTH))%SIMD_VECTOR_LENGTH;
    numXWithPad    = numx + padlength;   
    numVectorX = numXWithPad / SIMD_VECTOR_LENGTH;
    
    unsigned int numElementXYZ = numXWithPad * numy * numz;

    int maxnumthreads = 1;
    #ifdef _OPENMP
        maxnumthreads = omp_get_max_threads();
    #endif

    /* memory allocation */
    mask        = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    im_orig     = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    im_synth    = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*numElementXYZ,MEMORY_ALIGNMENT);
    //im_init_int = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //im_init_fz  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //im_end_int  = (SIMD_ELEMENT_TYPE *)_aligned_malloc(sizeof(SIMD_ELEMENT_TYPE)*SIMD_VECTOR_LENGTH*(numz+2),MEMORY_ALIGNMENT);
    //flag_invalid_allZ     =       (int *)_aligned_malloc(sizeof(int)*numVectorX*numy,MEMORY_ALIGNMENT);
    rss_current_mm     = (__MM *) _aligned_malloc(sizeof(__MM)*maxnumthreads,MEMORY_ALIGNMENT);
    rss_current_err_mm = (__MM *) _aligned_malloc(sizeof(__MM)*maxnumthreads,MEMORY_ALIGNMENT);


    /* SIMDベクトルの配列へ変換 */
    mask_mm     = (__MM *)mask;
    im_orig_mm  = (__MM *)im_orig;
    im_synth_mm  = (__MM *)im_synth;
    //imInitInt_mm = (__MM *)im_init_int;
    //imInitFz_mm  = (__MM *)im_init_fz;
    //imEndInt_mm  = (__MM *)im_end_int;
    

    ///* im_initを初期化 */
    //for (unsigned int c=0; c<numz+2; c++) {
    //    imInitInt_mm[c] = zeros_mm;
    //    imInitFz_mm[c]  = zeros_mm;
    //    imEndInt_mm[c]  = mask_true_mm;
    //}

    /* table_mask_xy を初期化 */
    table_mask_xy.resize(numy);
    for (unsigned int c=0; c<table_mask_xy.size(); c++) {
        table_mask_xy[c].clear();
        table_mask_xy[c].reserve(numVectorX/4);
    }
    
    /* table_elipsoid_xy を初期化 */
    table_elipsoid_xy.resize(numVectorX*numy);
    for (unsigned int c=0; c<table_elipsoid_xy.size(); c++) {
        table_elipsoid_xy[c].resize(1);
    }

    /* table_elipsoid_z を初期化 */
    table_elipsoid_z.resize(numz);
    for (unsigned int c=0; c<table_elipsoid_z.size(); c++) {
        table_elipsoid_z[c].resize(1);
    }

    /* 入力ポインタの実部ポインタの取得 */
	double *pim;
    pim = mxGetPr(prhs[1]); // im


    /* copy original image data to im_orig */        
    unsigned int cvxyz, cin;
    __MM tmpsqr_mm, mask_padding_mm;
    __MM thrint_mm = _MM_SET1( (SIMD_ELEMENT_TYPE) thrint );
    __MM rss_orig_mm = zeros_mm;
    __MM rss_orig_err_mm = zeros_mm;       
    
    // make mask for padding
    SIMD_ELEMENT_TYPE *mask_padding = (SIMD_ELEMENT_TYPE *)&mask_padding_mm;
    for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
        mask_padding[ce] = mask_true;
    }
    for (unsigned int ce=SIMD_VECTOR_LENGTH-padlength; ce<SIMD_VECTOR_LENGTH; ce++) {
        mask_padding[ce] = mask_false;
    }

    
    for (unsigned int cy=0; cy<numy; cy++) {        
        
        // for elements without padding        
        for (unsigned int cvx=0; cvx<numVectorX-1; cvx++) {
            __MM tmpmask;
            //__MM tmpmask_sum = mask_false_mm;
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cin   = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = (SIMD_ELEMENT_TYPE) pim[ cin + ce ];                          
                }
                tmpmask = _MM_CMP(im_orig_mm[cvxyz],thrint_mm,_CMP_GE_OS);
                im_orig_mm[cvxyz] = _MM_AND(im_orig_mm[cvxyz],tmpmask);
                //mask_mm[cvxyz] = tmpmask;
                mask_mm[cvxyz] = mask_true_mm;
                //tmpmask_sum = _MM_OR(tmpmask_sum,tmpmask);
                tmpsqr_mm = _MM_SQR(_MM_AND(tmpmask,im_orig_mm[cvxyz]));
                twosum(&rss_orig_mm,&rss_orig_err_mm,&tmpsqr_mm);
            }
            //flag_invalid_allZ[cy*numVectorX+cvx] = _MM_TESTZ1(tmpmask_sum);
            //if (!_MM_TESTZ1(tmpmask_sum)) {
                table_mask_xy[cy].push_back(cvx);
            //}
            
        }

        // for elements with padding
        unsigned int cvx = numVectorX - 1; {
            __MM tmpmask;
            //__MM tmpmask_sum = mask_false_mm;
            for (unsigned int cz=0; cz<numz; cz++) {
                cvxyz = numz*numVectorX*cy + numz*cvx + cz;
                cin   = numx*numy*cz + numx*cy + SIMD_VECTOR_LENGTH*cvx;
                for (unsigned int ce=0; ce<SIMD_VECTOR_LENGTH-padlength; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = (SIMD_ELEMENT_TYPE) pim[ cin + ce ];                          
                }
                for (unsigned int ce=SIMD_VECTOR_LENGTH-padlength; ce<SIMD_VECTOR_LENGTH; ce++) {
                    im_orig[ cvxyz*SIMD_VECTOR_LENGTH + ce ] = 0; // padding with 0
                }            
                tmpmask = _MM_AND(mask_padding_mm,_MM_CMP(im_orig_mm[cvxyz],thrint_mm,_CMP_GE_OS));
                im_orig_mm[cvxyz] = _MM_AND(im_orig_mm[cvxyz],tmpmask);
                //mask_mm[cvxyz] = tmpmask;
                mask_mm[cvxyz] = mask_padding_mm;
                //tmpmask_sum = _MM_OR(tmpmask_sum,tmpmask);
                tmpsqr_mm = _MM_SQR(_MM_AND(tmpmask,im_orig_mm[cvxyz]));
                twosum(&rss_orig_mm,&rss_orig_err_mm,&tmpsqr_mm);        
            }
            //flag_invalid_allZ[cy*numVectorX+cvx] = _MM_TESTZ1(tmpmask_sum);
            //if (!_MM_TESTZ1(tmpmask_sum)) {
                table_mask_xy[cy].push_back(cvx);
            //}
        }
    }

    SIMD_ELEMENT_TYPE *rss_orig     = (SIMD_ELEMENT_TYPE *)&rss_orig_mm;
    SIMD_ELEMENT_TYPE *rss_orig_err = (SIMD_ELEMENT_TYPE *)&rss_orig_err_mm;
    SIMD_ELEMENT_TYPE tmpsum = 0.0;
    SIMD_ELEMENT_TYPE tmperr = 0.0;
    for (unsigned int counter=0; counter<SIMD_VECTOR_LENGTH; counter++) {
        twosum(&tmpsum,&tmperr,&rss_orig[counter]);
        twosum(&tmpsum,&tmperr,&rss_orig_err[counter]);
    }
    tmpsum += tmperr;
    rss_orig_inv = 1.0/((double)tmpsum);    

}




void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

    #ifdef MY_ASSERT
    /* 変数宣言 */
    unsigned __int64 t_start, t_end;

    t_start = __rdtsc();
    #endif

	/* 入出力引数の数チェック */
	if (nrhs > 3 || nrhs < 1) {
		mexErrMsgTxt("入力引数の数は1から3までです。");
	} else if (nlhs > 3) {
		mexErrMsgTxt("出力引数が多すぎます。");
	}

	/* 入力引数が2つ以上のときは初期化を行ってから実際の計算を行う。
	 入力引数が1つのときは初期化をせずそのまま実際の計算を行う。 */
    if(nrhs > 1) { 
        initialize(nrhs,prhs); /* メモリの確保や元画像の格納などホスト側の初期化 */
    }

    //////#ifdef _OPENMP
    //////// set openmp stack size through KMP_STACKSIZE environment variable.
    ////////if (nlhs==3) {
    //////    int hoge = kmp_get_stacksize_s();
    //////    kmp_set_stacksize_s(numXWithPad*numy*numz*sizeof(SIMD_ELEMENT_TYPE)+16*1024*1024);
    //////    hoge = kmp_get_stacksize_s();
    ////////}
    //////#endif


    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("initialize: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif


    /* パラメータを格納 */
    setparam_init(prhs[0]);

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("setparam_init: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif

    /* EM-like optimization */
    em_optim();

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("em_optim: %I64d[clocks]\n", t_end-t_start);
    t_start = __rdtsc();
    #endif


    switch (nlhs) {
    case 3 :
        plhs[2] = mxCreateNumericArray(numdims,dims,mxDOUBLE_CLASS,mxREAL);
        setsynth(mxGetPr(plhs[2]));        
    case 2 :
        plhs[1] = mxCreateDoubleScalar(score);
    case 1 :        
        plhs[0] = mxCreateDoubleMatrix(NUM_PARAMS_INPUT,numelipsoids,mxREAL);
        setresult(mxGetPr(plhs[0]));        
    }

    #ifdef MY_ASSERT
    t_end = __rdtsc();
    mexPrintf("setresult: %I64d[clocks]\n", t_end-t_start);
    #endif
}


