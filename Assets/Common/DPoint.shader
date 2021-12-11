// Upgrade NOTE: replaced 'mul(UNITY_MATRIX_MVP,*)' with 'UnityObjectToClipPos(*)'

Shader "Instanced/DPoint" {
    Properties {
        _Colour ("Colour", COLOR) = (1, 1, 1, 1)
        _Size ("Size", float) = 0.5
        _Shininess("Shininess", Float) = 10 //Shininess
        _SpecColor("Specular Color", Color) = (1, 1, 1, 1) //Specular highlights color

    }

    SubShader {
        Pass {
            Tags { "LightMode"="ForwardBase" "Queue" = "Gemometry" "RenderType" = "Opaque" "IgnoreProjector" = "True" }
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
                float3 normal : NORMAL;
                float4 world_pos : TEXCOORD1;
                float2 uv : TEXCOORD0;
                float2 screen_pos : TEXCOORD2;
            };

            float _Size;
            float4 _Colour;

            StructuredBuffer<Particle> particle_buffer;
            sampler2D _CameraDepthTexture;
            //sampler2D  _CustomDepthTexture;

            v2f vert (appdata_full v, uint instanceID : SV_InstanceID) {
                // take in data from the compute buffer, filled with data each frame in SimRenderer
                // offsetting and scaling it from the (0...grid_res, 0...grid_res) resolution of our sim into a nicer range for rendering


                //float4 data = float4((particle_buffer[instanceID].x.xyz - float3(32, 32, 32)) * 0.1, 1.0);
                float4 data = float4(particle_buffer[instanceID].x.xyz, 1.0);
                // Scaling vertices by our base size param (configurable in the material) and the mass of the particle
                float3 localPosition = v.vertex.xyz * (_Size * data.w) * 10;
                float3 worldPosition = data.xyz + localPosition;

                // project into camera space
                v2f o;
                o.world_pos = float4(worldPosition, 1.0);
                o.pos = UnityObjectToClipPos(worldPosition);
                o.uv = v.texcoord;
                o.normal = normalize(mul(float4(v.normal, 0.0), unity_WorldToObject).xyz); //Calculate the normal
                o.screen_pos = ComputeScreenPos(o.pos);
                //o.screenPos = ComputeScreenPos(o.pos);
                return o;
            }

            fixed4 frag(v2f i) : SV_Target{
                //float d = length(_WorldSpaceCameraPos - i.world_pos);
                float d = Linear01Depth(tex2D(_CameraDepthTexture, i.screen_pos.xy).r);
                //float d = tex2D(_CustomDepthTexture, i.screen_pos.xy).r;
                //float3 n;
                //float d;
                //DecodeDepthNormal(tex2D(_CameraDepthNormalsTexture, i.screen_pos.xy),d,n);
                
                //float4 col = float4(d / 50, d / 50, d / 50 + .5, 1);
                //float4 col = float4(.4,.4,.9, .7);
                float4 col = float4(d,d,d, 1);
                //return float4(d/50,d/50,d/50, 1);
                float3 normalDirection = normalize(i.normal);
                float3 viewDirection = normalize(_WorldSpaceCameraPos - i.world_pos.xyz);

                float3 vert2LightSource = _WorldSpaceLightPos0.xyz - i.world_pos.xyz;
                float oneOverDistance = 1.0 / length(vert2LightSource);
                float attenuation = lerp(1.0, oneOverDistance, _WorldSpaceLightPos0.w); //Optimization for spot lights. This isn't needed if you're just getting started.
                float3 lightDirection = _WorldSpaceLightPos0.xyz - i.world_pos.xyz * _WorldSpaceLightPos0.w;

                float3 ambientLighting = UNITY_LIGHTMODEL_AMBIENT.rgb * col.rgb; //Ambient component
                float3 diffuseReflection = attenuation * _LightColor0.rgb * col.rgb * max(0.0, dot(normalDirection, lightDirection)); //Diffuse component
                float3 specularReflection;
                if (dot(i.normal, lightDirection) < 0.0) //Light on the wrong side - no specular
                {
                    specularReflection = float3(0.0, 0.0, 0.0);
                }
                else
                {
                    //Specular component
                    specularReflection = attenuation * _LightColor0.rgb * _SpecColor.rgb * pow(max(0.0, dot(reflect(-lightDirection, normalDirection), viewDirection)), _Shininess);
                }
                float3 color = (ambientLighting + diffuseReflection) * col +specularReflection; //Texture is not applient on specularReflection
                return float4(color, 1.0);
            }
            ENDCG
        }
        Pass
        {
            Tags{ "LightMode" = "ShadowCaster" }
            Name "Shadow Cast"

            CGPROGRAM
            #pragma vertex vert
            #pragma fragment frag
            #pragma multi_compile_shadowcaster

            #include "UnityCG.cginc"

            struct v2f
            {
                V2F_SHADOW_CASTER;
            };

            v2f vert(appdata_base v)
            {
                v2f o;
                TRANSFER_SHADOW_CASTER_NORMALOFFSET(o);
                return o;
            }

            fixed4 frag(v2f i) : SV_Target
            {
                SHADOW_CASTER_FRAGMENT(i);
            }
            ENDCG
        }
    }
}