// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Instanced/DepthShader" {
    Properties{
        _Colour("Colour", COLOR) = (1, 1, 1, 1)
        _Size("Size", float) = 0.5
        _Shininess("Shininess", Float) = 10 //Shininess
        _SpecColor("Specular Color", Color) = (1, 1, 1, 1) //Specular highlights color

    }

        SubShader{
            Tags {"RenderType" = "Opaque"}
            Pass {
                ZWrite On
                Blend SrcAlpha OneMinusSrcAlpha

                CGPROGRAM

                #pragma glsl
                #pragma vertex vert
                #pragma fragment frag
                #pragma multi_compile_fwdbase nolightmap nodirlightmap nodynlightmap novertexlight
                #pragma target 4.5

                #include "UnityCG.cginc"
                uniform float4 _LightColor0;
                uniform float4 _SpecColor;
                uniform float _Shininess;
                // matches the structure of our data on the CPU side
                struct Particle {
                    float3 x;
                    float3 v;
                    float3x3 C;
                    float mass;
                    float padding;
                };

                /*struct v2f {
                    float4 pos : SV_POSITION;
                    float4 world_pos: TEXCOORD0;
                    //float4 screenPos : TEXCOORD0;
                };*/

                struct v2f {
                    float4 pos : SV_POSITION;
                    float4 world_pos : TEXCOORD1;
                };

                float _Size;
                float4 _Colour;

                StructuredBuffer<Particle> particle_buffer;
                sampler2D _CameraDepthTexture;

                v2f vert(appdata_full v, uint instanceID : SV_InstanceID) {
                    // take in data from the compute buffer, filled with data each frame in SimRenderer
                    // offsetting and scaling it from the (0...grid_res, 0...grid_res) resolution of our sim into a nicer range for rendering


                    //float4 data = float4((particle_buffer[instanceID].x.xyz - float3(32, 32, 32)) * 0.1, 1.0);
                    float4 data = float4(particle_buffer[instanceID].x.xyz, 1.0);
                    // Scaling vertices by our base size param (configurable in the material) and the mass of the particle
                    float3 localPosition = v.vertex.xyz * (_Size * data.w) * 10;
                    float3 worldPosition = data.xyz + localPosition;
                    o.pos = UnityObjectToClipPos(worldPosition);
                    // project into camera space
                    v2f o;
                    o.world_pos = float4(worldPosition, 1.0);
                    return o;
                }

                fixed4 frag(v2f i) : SV_Target{
                    float d = length(_WorldSpaceCameraPos - i.world_pos);
                    //float d = Linear01Depth(tex2D(_CameraDepthTexture, i.screen_pos.xy).r);
                    return float4(d,1,d, 1.0);
                }
            ENDCG
            }
    }
}